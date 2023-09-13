#include <hibf/config.hpp>

int main()
{
    auto my_input = [&](size_t const /* user_bin_id */, seqan::hibf::insert_iterator it) // fixed parameters!
    {
        it = 42; // assign something that is convertible to uint64_t
    };

    seqan::hibf::config config{.input_fn = my_input, .number_of_user_bins = 12};
}
