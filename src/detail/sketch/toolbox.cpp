#include <cassert>    // for assert
#include <cmath>      // for floor, sqrt
#include <cstddef>    // for size_t
#include <functional> // for greater
#include <limits>     // for numeric_limits
#include <omp.h>      // for omp_get_thread_num
#include <utility>    // for move
#include <vector>     // for vector

#include <hibf/contrib/robin_hood.hpp>        // for unordered_flat_map, pair
#include <hibf/detail/sketch/hyperloglog.hpp> // for hyperloglog
#include <hibf/detail/sketch/toolbox.hpp>     // for toolbox

namespace hibf::sketch
{

void toolbox::cluster_bins(std::vector<size_t> & permutation,
                           size_t const first,
                           size_t const last,
                           size_t const num_threads)
{
    assert(num_threads >= 1);
    assert(positions->size() <= sketches->size());
    assert((first == 0) == permutation.empty());

    size_t const n = sketches->size();
    size_t const chunk_size = std::floor(std::sqrt(n));

    size_t const prune_steps = chunk_size;
    size_t steps_without_prune = 0;

    size_t const none = std::numeric_limits<size_t>::max();
    /* internal map that stores the distances
        *
        * The first layer is a hash map with the ids of active clusters as keys.
        * The values (second layer) are priority queues with neighbors of the cluster
        * with the respective key in the first layer.
        * These neighbors are themselves clusters with an id and store a distance to the
        * cluster of the first layer.
        */
    distance_matrix dist;
    dist.reserve(n + 1);

    // map that indicates which ids of clusters are still in the distance matrix
    // the values are the indices where the priority queue for the given id as key can be found in dist
    robin_hood::unordered_flat_map<size_t, size_t> remaining_ids;

    // clustering tree stored implicitly in a vector
    std::vector<clustering_node> clustering;
    clustering.reserve(2 * n);

    // cache for hll cardinality estimates
    std::vector<double> estimates;
    estimates.reserve(2 * n);

    // every thread will write its observed id with minimal distance to some other here
    // id == none means that the thread observed only empty or no priority queues
    std::vector<size_t> min_ids(num_threads, none);

    // these will be the new ids for new clusters
    // the first one is invalid, but it will be incremented before it is used for the first time
    size_t new_id = last - 1;

    // initialize clustering and estimates
    for (size_t id = first; id < last; ++id)
    {
        // id i is at the index i - first
        clustering.push_back({none, none, (*sketches)[(*positions)[id]]});
        estimates.emplace_back((*sketches)[(*positions)[id]].estimate());
    }

    // if this is not the first group, we want to have one overlapping bin
    size_t previous_rightmost = none;
    if (first != 0)
    {
        // For all other clusters, their id is also their index in filesnames, sketches etc. .
        // This is important, because their id is then inserted into the clustering.
        // This does not work for previous rightmost, because its index does not necessarily lie on
        // the continuous spectrum from first to last. We run into a problem, because the entries are
        // stored in vectors. Therefore we give previous_rightmost a different id (==last). This is
        // fine, because we only need the HLL sketch of the actual index. previous_rightmost will be ignored
        // in the traceback anyway and won't be added to the permutation in this step.
        size_t const actual_previous_rightmost = permutation.back();
        ++new_id;
        previous_rightmost = new_id;

        clustering.push_back({none, none, (*sketches)[(*positions)[actual_previous_rightmost]]});
        estimates.emplace_back((*sketches)[(*positions)[actual_previous_rightmost]].estimate());
    }

    // initialize priority queues in the distance matrix (sequentially)
    for (size_t id = first; id < first + clustering.size(); ++id)
    {
        // empty priority queue for every item in clustering
        dist.push_back({id, prio_queue{}});
        remaining_ids[id] = id - first;
    }

#pragma omp parallel num_threads(num_threads)
    {
        double min_dist = std::numeric_limits<double>::max();
// minimum distance exclusively for this thread

// initialize all the priority queues of the distance matrix
// while doing that, compute the first min_id
#pragma omp for schedule(nonmonotonic : dynamic, chunk_size)
        for (size_t i = 0; i < clustering.size(); ++i)
        {
            for (size_t j = 0; j < clustering.size(); ++j)
            {
                // we only want one diagonal of the distance matrix
                if (i < j)
                {
                    // this must be a copy, because merging changes the hll sketch
                    hyperloglog temp_hll = clustering[i].hll;
                    double const estimate_ij = temp_hll.merge_and_estimate_SIMD(clustering[j].hll);
                    // Jaccard distance estimate
                    double const distance = 2 - (estimates[i] + estimates[j]) / estimate_ij;
                    dist[i].pq.push({j + first, distance});
                }
            }
            if (dist[i].pq.empty())
                continue;

            // check if the just initialized priority queue contains the minimum value for this thread
            neighbor const & curr = dist[i].pq.top();
            if (curr.dist < min_dist)
            {
                min_dist = curr.dist;
                min_ids[omp_get_thread_num()] = dist[i].id;
            }
        } // implicit barrier

// a single thread shuffles dist to approximately balance loads in static scheduling
#pragma omp single
        random_shuffle(dist, remaining_ids);

        // main loop of the clustering
        // keep merging nodes until we have a complete tree
        while (remaining_ids.size() > 1)
        {
// Wait for all threads to have evaluated remaining_ids.size() as remaining_ids
// may be modified by the following pragma omp single.
#pragma omp barrier

#pragma omp single
            {
                // perform critical update
                // increment id for the new cluster (must be done at the beginning)
                ++new_id;

                // compute the final min_id from the min_ids of the worker threads
                size_t min_id = min_ids[0];
                double min_dist = std::numeric_limits<double>::max();
                for (auto candidate_id : min_ids)
                {
                    // check if the thread saw any id
                    if (candidate_id == none)
                        continue;

                    size_t const dist_index = remaining_ids.at(candidate_id);
                    neighbor const & curr = dist[dist_index].pq.top();
                    if (curr.dist < min_dist)
                    {
                        min_dist = curr.dist;
                        min_id = candidate_id;
                    }
                }

                size_t const min_index = remaining_ids.at(min_id); // how can min_id be none?
                size_t const neighbor_id = dist[min_index].pq.top().id;

                // merge the two nodes with minimal distance together insert the new node into the clustering
                clustering.push_back({min_id, neighbor_id, std::move(clustering[min_id - first].hll)});
                estimates.emplace_back(
                    clustering.back().hll.merge_and_estimate_SIMD(clustering[neighbor_id - first].hll));

                // remove old ids
                remaining_ids.erase(min_id);
                remaining_ids.erase(neighbor_id);

                // overwrite one of the removed entries with the new one
                remaining_ids[new_id] = min_index;
                dist[min_index] = {new_id, prio_queue{}};

                // prune the distance matrix to reduce overhead due to inactive entries
                ++steps_without_prune;
                if (steps_without_prune > prune_steps)
                {
                    prune(dist, remaining_ids);
                    steps_without_prune = 0;
                }
            } // implicit barrier

            // reset values for the computation of the new minimum
            min_ids[omp_get_thread_num()] = none;
            min_dist = std::numeric_limits<double>::max();

            hyperloglog const new_hll = clustering.back().hll;

// update distances in dist
// while doing that, compute the new min_id
#pragma omp for schedule(static)
            for (size_t i = 0; i < dist.size(); ++i)
            {
                size_t other_id = dist[i].id;
                if (other_id == new_id || !remaining_ids.contains(other_id))
                    continue;

                // this must be a copy, because merge_and_estimate_SIMD() changes the hll
                hyperloglog temp_hll = new_hll;
                double const estimate_ij = temp_hll.merge_and_estimate_SIMD(clustering[other_id - first].hll);
                // Jaccard distance estimate
                double const distance = 2 - (estimates[other_id - first] + estimates.back()) / estimate_ij;
                dist[i].pq.push({new_id, distance});

                // make sure the closest neighbor is not yet deleted (this is a lazy update)
                while (!remaining_ids.contains(dist[i].pq.top().id))
                {
                    dist[i].pq.pop();
                }

                // check if the just updated priority queue contains the minimum value for this thread
                neighbor const & curr = dist[i].pq.top();
                if (curr.dist < min_dist)
                {
                    min_dist = curr.dist;
                    min_ids[omp_get_thread_num()] = other_id;
                }
            } // implicit barrier
        }
    } // end of the parallel region

    size_t final_root_index = remaining_ids.begin()->second;
    size_t final_root_id = dist[final_root_index].id;

    // rotate the previous rightmost to the left so that it has the correct place in the permutation
    if (first != 0)
    {
        rotate(clustering, previous_rightmost, first, final_root_id);
    }

    // traceback into permutation and ignore the previous rightmost
    trace(clustering, permutation, previous_rightmost, first, final_root_id);
}

} // namespace hibf::sketch
