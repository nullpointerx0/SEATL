#ifndef SEATL_SOLVER_HPP
#define SEATL_SOLVER_HPP

/*
 * SEATL Solver - Super (a,d)-Edge-Antimagic Total Labeling
 *
 * For complete bipartite graph K_{m,n}:
 *   - Vertex labels: {1, 2, ..., m+n}, each used once
 *   - Edge labels:   {m+n+1, ..., m+n+mn}, each used once
 *   - Edge weight:   w(e) = label(u) + label(v) + label(e)
 *   - Goal: all edge weights form arithmetic progression {a, a+1, ..., a+mn-1}
 *
 * 7-Layer Architecture:
 *   Layer 1: Master formula (algebraic pruning)
 *   Layer 2: AC-3 constraint propagation
 *   Layer 3: Greedy initialization
 *   Layer 4: Iterative deepening DFS
 *   Layer 5: Most-constrained-variable ordering
 *   Layer 6: Memoization
 *   Layer 7: Restart with diversification
 */

#include <vector>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <optional>
#include <chrono>
#include <string>
#include <utility>

namespace seatl {

// ============================================================
// Configuration
// ============================================================

struct Config {
    // Layer toggles
    bool use_master_formula      = true;   // Layer 1
    bool use_ac3                 = true;   // Layer 2
    bool use_greedy              = true;   // Layer 3
    bool use_iterative_deepening = true;   // Layer 4
    bool use_mcv_ordering        = true;   // Layer 5
    bool use_memoization         = true;   // Layer 6
    bool use_restart             = true;   // Layer 7

    // Extra strategies (Layers 8-11)
    bool use_lcv_value_ordering  = true;   // L8: try least-constraining label first
    bool use_target_propagation  = true;   // L9: prune based on target slot uniqueness
    bool use_partition_ordering  = true;   // L10: try most-promising partitions first
    bool use_local_search        = true;   // L11: min-conflicts repair as fallback

    // Local search tuning
    int    local_search_iterations = 10000;  // Max repair steps per attempt
    int    local_search_attempts   = 3;      // Random restarts within local search

    // Enumeration mode: find all solutions instead of stopping at first.
    // When true, the solver runs DFS in enumeration mode (local search is skipped
    // because it can't enumerate completely). For large graphs this is slow; use
    // max_solutions and time_limit_seconds to bound the work.
    bool enumerate_all       = false;
    int  max_solutions       = 100;        // Hard cap on solutions returned
    bool dedupe_by_canonical = true;       // Skip M<->N mirror duplicates

    // Tuning
    int    max_dfs_depth        = 50;       // Hard recursion ceiling
    int    id_initial_depth     = 10;       // First iterative-deepening cap
    int    id_depth_step        = 10;       // How much each retry adds
    int    restart_attempts     = 5;        // Number of restart tries
    int    memoization_limit    = 100000;   // Max memo table entries
    double time_limit_seconds   = 15.0;     // Wall-clock budget (total)
    double per_partition_time   = 1.0;      // Wall-clock budget (per partition)

    // Diagnostics
    bool verbose      = false;
    bool stop_on_first = true;  // Return on first solution found
};

// ============================================================
// Result
// ============================================================

struct Result {
    bool found      = false;
    bool timed_out  = false;
    int  a          = -1;       // Starting weight of AP

    std::vector<int> set_m;     // Vertex labels in M
    std::vector<int> set_n;     // Vertex labels in N

    // Each entry: ((vertex_in_M, vertex_in_N), edge_label)
    std::vector<std::pair<std::pair<int,int>, int>> edge_labeling;
    std::vector<int> edge_weights;  // Sorted; should be {a, a+1, ..., a+mn-1}

    // Diagnostics
    long long partitions_tested = 0;
    long long partitions_pruned_layer1 = 0;
    long long partitions_pruned_layer2 = 0;
    long long greedy_successes  = 0;
    long long dfs_invocations   = 0;
    long long memo_hits         = 0;
    long long restarts_used     = 0;
    double    elapsed_ms        = 0.0;
    std::string status;
};

// A single SEATL solution (slimmer than Result; used for enumeration).
struct Solution {
    int a = -1;
    std::vector<int> set_m;
    std::vector<int> set_n;
    std::vector<std::pair<std::pair<int,int>, int>> edge_labeling;
    std::vector<int> edge_weights;
};

// Container returned by enumerate_solutions().
struct EnumerationResult {
    std::vector<Solution> solutions;
    bool        complete      = false;  // True iff every partition was exhausted
    bool        timed_out     = false;
    bool        capped        = false;  // True iff max_solutions was hit
    long long   partitions_tested = 0;
    long long   partitions_pruned = 0;
    double      elapsed_ms    = 0.0;
    std::string status;
};

// ============================================================
// Solver
// ============================================================

class Solver {
public:
    Solver(int m, int n, const Config& cfg = Config());

    // Main entry: try every canonical partition, return first solution
    Result solve();

    // Enumeration entry: find all SEATL solutions (up to max_solutions / time limit).
    // Uses systematic DFS only (no local search), enumerating every leaf.
    EnumerationResult solve_all();

private:
    // ---- Problem parameters ----
    const int m_, n_, mn_;
    const int min_label_;       // m + n + 1   (smallest edge label)
    const int max_label_;       // m + n + mn  (largest edge label)
    const Config cfg_;

    // ---- Timing ----
    std::chrono::steady_clock::time_point start_time_;
    std::chrono::steady_clock::time_point partition_deadline_;
    bool has_partition_deadline_ = false;
    bool time_exceeded() const;       // true if total OR per-partition deadline hit
    bool total_time_exceeded() const; // true only if total budget exhausted

    // ---- Layer 1: Master Formula ----
    // Returns 'a' if integer, -1 otherwise.
    int compute_starting_weight(const std::vector<int>& M,
                                const std::vector<int>& N) const;

    // Sum of all edge labels {min_label_..max_label_}
    long long edge_label_sum() const;

    // ---- Layer 2: AC-3 ----
    struct ConstraintState {
        // feasible[i] = set of edge-labels that could be assigned to edge i
        std::vector<std::set<int>> feasible;
        bool consistent = true;
    };

    // Build initial feasibility table; returns false if any edge has empty domain.
    bool build_initial_feasibility(const std::vector<int>& pair_sums,
                                   const std::vector<int>& targets,
                                   ConstraintState& out) const;

    // Run AC-3 fixed-point loop. Returns false if inconsistent.
    bool ac3_propagate(const std::vector<int>& pair_sums,
                       const std::vector<int>& targets,
                       ConstraintState& state) const;

    // Forward checking during DFS: remove `assigned_label` from all unassigned edges.
    // Records changes for undo. Returns false if any unassigned edge becomes empty.
    bool forward_check(int assigned_label,
                       const std::vector<bool>& edge_assigned,
                       ConstraintState& state,
                       std::vector<std::pair<int,int>>& removals) const;

    void undo_forward_check(ConstraintState& state,
                            const std::vector<std::pair<int,int>>& removals) const;

    // ---- Layer 3: Greedy ----
    // Tries simple in-order assignment; returns assignment vector or nullopt.
    std::optional<std::vector<int>> greedy_assign(
        const std::vector<int>& pair_sums,
        const ConstraintState& initial_state) const;

    // ---- Layers 4 & 5: Iterative Deepening + MCV ----
    bool iterative_deepening_search(
        const std::vector<int>& pair_sums,
        ConstraintState& state,
        std::vector<int>& assignment,
        const std::vector<int>& edge_order_hint);

    // Returns: 0 = success, 1 = failure, 2 = depth-limit hit
    enum class DfsOutcome { Success, Failure, DepthExceeded };

    DfsOutcome dfs(int depth,
                   int depth_limit,
                   const std::vector<int>& pair_sums,
                   ConstraintState& state,
                   std::vector<int>& assignment,
                   std::vector<bool>& edge_assigned,
                   const std::vector<int>* fixed_order);

    // Enumeration variant: records every complete assignment to `solutions_out`
    // instead of returning on the first one. Stops once max_solutions is reached
    // or time budget is exhausted. Returns true iff search was fully exhaustive.
    bool dfs_enumerate(int depth,
                       int depth_limit,
                       const std::vector<int>& pair_sums,
                       const std::vector<int>& M,
                       const std::vector<int>& N,
                       int a,
                       ConstraintState& state,
                       std::vector<int>& assignment,
                       std::vector<bool>& edge_assigned,
                       std::vector<Solution>& solutions_out,
                       int max_solutions);

    // Layer 5: pick next edge to assign (most constrained = fewest options).
    int select_next_edge_mcv(const ConstraintState& state,
                             const std::vector<bool>& edge_assigned) const;

    // ---- Layer 6: Memoization ----
    // Key: pair (edge_idx, hash of used-labels bitset). Value: known-failure flag.
    std::unordered_set<unsigned long long> memo_failures_;
    long long memo_hits_count_ = 0;

    unsigned long long memo_key(int edge_idx,
                                const std::vector<int>& assignment) const;

    bool memo_check_failure(int edge_idx,
                            const std::vector<int>& assignment) const;

    void memo_record_failure(int edge_idx,
                             const std::vector<int>& assignment);

    // ---- Layer 7: Restart ----
    bool solve_with_restarts(const std::vector<int>& pair_sums,
                             ConstraintState& state,
                             std::vector<int>& assignment,
                             int& restarts_used);

    // ---- Layer 11: Local search (min-conflicts) ----
    // Independent of DFS; repairs random initial assignments toward a valid SEATL.
    bool local_search(const std::vector<int>& pair_sums,
                      const std::vector<int>& targets,
                      const ConstraintState& state,
                      std::vector<int>& assignment_out);

    // ---- Helpers ----
    // Compute pair sums for K_{m,n} given partition (M, N).
    // Returns vector of size mn; element i corresponds to (M[i/n], N[i%n]).
    std::vector<int> compute_pair_sums(const std::vector<int>& M,
                                       const std::vector<int>& N) const;

    // Targets are {a, a+1, ..., a+mn-1}.
    std::vector<int> compute_targets(int a) const;

    // Generate all canonical partitions (M contains label 1 to break symmetry).
    void generate_partitions(
        std::vector<std::pair<std::vector<int>, std::vector<int>>>& out) const;

    // Build the final Result struct from a successful assignment.
    Result build_result(const std::vector<int>& M,
                        const std::vector<int>& N,
                        const std::vector<int>& assignment,
                        int a) const;
};

} // namespace seatl

#endif // SEATL_SOLVER_HPP
