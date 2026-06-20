#ifndef SEATL_SOLVER_XL_HPP
#define SEATL_SOLVER_XL_HPP

/*
 * SEATL Solver XL — Scales to K_{20,20} in minutes (not hours).
 *
 * 12-Layer Architecture (extends 11-layer):
 *   Layers 1-11: same as before (master formula, AC-3, greedy, DFS, MCV,
 *                memoization, restart, LCV, target propagation, partition
 *                ordering, local search)
 *   Layer 12:    Parallel Tabu Search with bitset assignments
 *                + smart partition seeding (skip enumeration for large m,n)
 *
 * For large graphs (m+n >= 20), the partition enumeration of ~10^11 candidates
 * is the bottleneck, not the per-partition search. The XL solver bypasses this
 * by seeding directly from algebraically promising partitions:
 *
 *   1. SKIP partition enumeration when |M|*|N| > 200
 *   2. Start from the "balanced block" partition M={1..m}, N={m+1..m+n}
 *   3. Run parallel tabu search on multiple threads (OpenMP)
 *   4. If primary partition fails, jitter to a neighbor partition
 *
 * Performance target: K_{20,20} in < 5 minutes on 8-core CPU.
 */

#include "seatl_solver.hpp"

namespace seatl {

// XL-specific configuration
struct ConfigXL : public Config {
    // Layer 12: Parallel Tabu Search
    bool use_parallel_tabu     = true;     // Enable OpenMP parallelism
    int  num_threads           = 0;        // 0 = auto-detect (use all cores)
    int  tabu_iterations       = 500000;   // Max iterations per attempt
    int  tabu_tenure           = 50;       // How long a move is forbidden
    int  tabu_attempts         = 16;       // Parallel restarts

    // Smart partition seeding
    bool skip_enumeration_xl   = true;     // For large m,n, don't enumerate
    int  enumeration_threshold = 30;       // m*n > this triggers XL mode (almost everything)
    int  partition_jitter      = 8;        // Try N neighboring partitions if needed

    // Aggressive optimization for large graphs
    bool use_bitset_repr       = true;     // Use bitsets instead of std::set
    int  partition_search_depth = 10;      // How many partition swaps to try
};

// XL Solver — drop-in replacement for Solver, optimized for large graphs
class SolverXL {
public:
    SolverXL(int m, int n, const ConfigXL& cfg = ConfigXL());

    // Main entry — automatically chooses strategy based on graph size
    Result solve();

private:
    const int m_, n_, mn_;
    const int min_label_, max_label_;
    const ConfigXL cfg_;

    std::chrono::steady_clock::time_point start_time_;

    bool time_exceeded() const;

    // Compute 'a' (master formula)
    int compute_a(const std::vector<int>& M, const std::vector<int>& N) const;

    // Compute pair sums for a partition
    std::vector<int> compute_pair_sums(const std::vector<int>& M,
                                       const std::vector<int>& N) const;

    // Layer 12: parallel tabu search on a single partition
    // Returns assignment (edge_idx -> label) if successful
    bool tabu_search(const std::vector<int>& pair_sums,
                     int a,
                     std::vector<int>& assignment_out,
                     int thread_seed);

    // Tabu search worker (called in parallel)
    bool tabu_worker(const std::vector<int>& pair_sums,
                     int a,
                     std::vector<int>& best_assignment,
                     int seed,
                     std::chrono::steady_clock::time_point deadline);

    // Smart partition generation for XL mode
    void generate_seed_partitions(
        std::vector<std::pair<std::vector<int>, std::vector<int>>>& out) const;

    // For small/medium graphs, fall back to standard solver
    Result solve_with_standard_solver();

    // Build final result
    Result build_result(const std::vector<int>& M,
                        const std::vector<int>& N,
                        const std::vector<int>& assignment,
                        int a) const;
};

} // namespace seatl

#endif // SEATL_SOLVER_XL_HPP
