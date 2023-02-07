#!/usr/bin/env bash
set -Eeuo pipefail

args=${@/--branch-counts/""}
args=${args/--branch-probabilities/""}

exec gcov $args
