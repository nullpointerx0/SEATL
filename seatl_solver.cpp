/*
 * SEATL Solver implementation - all 7 layers.
 *
 * See seatl_solver.hpp for the architecture overview.
 */

#include "seatl_solver.hpp"

#include <algorithm>
#include <numeric>
#include <iostream>
#include <random>
#include <functional>
#include <cassert>
#include <climits>
#include <cstdlib>

namespace seatl {

// ============================================================
// Construction & timing
// ============================================================

Solver::Solver(int m, int n, const Config& cfg)
    : m_(m),
      n_(n),
      mn_(m * n),
      min_label_(m + n + 1),
      max_label_(m + n + m * n),
      cfg_(cfg) {
    // Defensive checks
    if (m < 2 || n < 2) {
        throw std::invalid_argument("SEATL requires m,n >= 2");
    }
}

bool Solver::time_exceeded() const {
    auto now = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(now - start_time_).count();
    if (elapsed > cfg_.time_limit_seconds) return true;
    if (has_partition_deadline_ && now > partition_deadline_) return true;
    return false;
}

bool Solver::total_time_exceeded() const {
    auto now = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(now - start_time_).count();
    return elapsed > cfg_.time_limit_seconds;
}

// ============================================================
// LAYER 1: Master Formula
// ============================================================
//
// Sum of all edge weights must equal sum of AP {a, a+1, ..., a+mn-1}
//   = mn*a + mn(mn-1)/2
//
// Sum of all edge weights, expanded by counting how often each label appears:
//   Each vertex in M appears in n edges  -> contributes n * sum(M)
//   Each vertex in N appears in m edges  -> contributes m * sum(N)
//   Each edge label appears once         -> contributes sum(edge labels)
//
// Setting equal and solving:
//   a = (n*sum_M + m*sum_N + sum_E - mn(mn-1)/2) / mn
//
// If 'a' is not an integer, the partition is provably infeasible.

long long Solver::edge_label_sum() const {
    // Sum of arithmetic series from min_label_ to max_label_.
    long long first = min_label_;
    long long last  = max_label_;
    long long count = mn_;
    return (first + last) * count / 2;
}

int Solver::compute_starting_weight(const std::vector<int>& M,
                                    const std::vector<int>& N) const {
    long long sum_m = std::accumulate(M.begin(), M.end(), 0LL);
    long long sum_n = std::accumulate(N.begin(), N.end(), 0LL);
    long long sum_e = edge_label_sum();

    long long numerator = (long long)n_ * sum_m
                        + (long long)m_ * sum_n
                        + sum_e
                        - (long long)mn_ * (mn_ - 1) / 2;

    if (numerator % mn_ != 0) return -1;
    long long a = numerator / mn_;
    if (a < 0) return -1;
    return (int)a;
}

// ============================================================
// LAYER 2: AC-3 Constraint Propagation
// ============================================================

bool Solver::build_initial_feasibility(const std::vector<int>& pair_sums,
                                       const std::vector<int>& targets,
                                       ConstraintState& out) const {
    out.feasible.assign(pair_sums.size(), {});
    out.consistent = true;

    // For edge i with pair sum P_i, target weight w requires edge label (w - P_i).
    // Label must be in [min_label_, max_label_].
    for (size_t i = 0; i < pair_sums.size(); ++i) {
        int p = pair_sums[i];
        for (int w : targets) {
            int needed = w - p;
            if (needed >= min_label_ && needed <= max_label_) {
                out.feasible[i].insert(needed);
            }
        }
        if (out.feasible[i].empty()) {
            out.consistent = false;
            return false;
        }
    }
    return true;
}

bool Solver::ac3_propagate(const std::vector<int>& pair_sums,
                           const std::vector<int>& targets,
                           ConstraintState& state) const {
    // Two propagation rules to fixed point:
    //   (R1) If edge i has only one feasible label L, then for every other edge j:
    //          - L cannot appear in j's feasible set (label uniqueness).
    //          - Any label in j's feasible set that produces the same target weight
    //            (P_i + L) must be removed (target weight uniqueness).
    //   (R2) For each target weight W, count how many (edge, label) pairs in the
    //          current feasibility table can produce W. If only one such pair exists,
    //          that edge is forced to that label.
    //
    // R2 is the "from the target side" rule that helps a lot when symmetry stops L1
    // from pruning anything by itself.

    bool changed = true;
    int iterations = 0;
    const int MAX_ITER = 50;

    while (changed && iterations < MAX_ITER) {
        changed = false;
        iterations++;
        if (time_exceeded()) return false;

        // ---- R1: singleton domains force removals elsewhere ----
        for (size_t i = 0; i < state.feasible.size(); ++i) {
            if (state.feasible[i].size() != 1) continue;
            int forced_label = *state.feasible[i].begin();
            int forced_weight = pair_sums[i] + forced_label;

            for (size_t j = 0; j < state.feasible.size(); ++j) {
                if (j == i) continue;
                // Remove forced_label
                {
                    auto it = state.feasible[j].find(forced_label);
                    if (it != state.feasible[j].end()) {
                        state.feasible[j].erase(it);
                        changed = true;
                        if (state.feasible[j].empty()) {
                            state.consistent = false;
                            return false;
                        }
                    }
                }
                // Remove labels in j's domain that would produce forced_weight
                int forbidden = forced_weight - pair_sums[j];
                if (forbidden != forced_label) {
                    auto it2 = state.feasible[j].find(forbidden);
                    if (it2 != state.feasible[j].end()) {
                        state.feasible[j].erase(it2);
                        changed = true;
                        if (state.feasible[j].empty()) {
                            state.consistent = false;
                            return false;
                        }
                    }
                }
            }
        }

        // ---- R2: each target must be claimable by exactly one edge ----
        // For each target W, scan: which (edge, label) pairs produce W?
        // If count is 0 -> infeasible. If count is 1 -> that edge forced.
        for (int w : targets) {
            int producer_edge = -1;
            int producer_label = -1;
            int count = 0;
            for (size_t i = 0; i < state.feasible.size(); ++i) {
                int needed = w - pair_sums[i];
                if (state.feasible[i].count(needed)) {
                    count++;
                    producer_edge = (int)i;
                    producer_label = needed;
                    if (count > 1) break;  // Can't force; move on
                }
            }
            if (count == 0) {
                state.consistent = false;
                return false;
            }
            if (count == 1) {
                // Force producer_edge to take producer_label.
                if (state.feasible[producer_edge].size() > 1) {
                    state.feasible[producer_edge].clear();
                    state.feasible[producer_edge].insert(producer_label);
                    changed = true;
                }
            }
        }
    }

    return state.consistent;
}

bool Solver::forward_check(int assigned_label,
                           const std::vector<bool>& edge_assigned,
                           ConstraintState& state,
                           std::vector<std::pair<int,int>>& removals) const {
    // We must enforce TWO uniqueness constraints when an edge is assigned:
    //   (a) The label itself can't be reused -> remove it from other domains.
    //   (b) The weight produced (pair_sum + label) is a target slot now claimed,
    //       so no other edge may produce that same weight.
    //
    // Constraint (a) is straightforward: erase `assigned_label` from each unassigned
    // edge's feasible set.
    //
    // Constraint (b) requires knowing which edge was just assigned (to compute the
    // weight). The caller passes `assigned_label` only - so we infer the produced
    // weight via the singleton domain we just set on edge `assigned_edge_idx`.
    // To keep the interface simple, we compute the produced weight by finding
    // which assigned edge has `assigned_label` in its (now singleton) domain
    // and using its pair_sum. This is wasteful; instead we rely on the caller
    // to track pair_sums separately. We'll handle that via a member call.
    //
    // Implementation: pass pair_sums via a member pointer set before each DFS.
    for (size_t i = 0; i < state.feasible.size(); ++i) {
        if (edge_assigned[i]) continue;
        auto it = state.feasible[i].find(assigned_label);
        if (it != state.feasible[i].end()) {
            state.feasible[i].erase(it);
            removals.push_back({(int)i, assigned_label});
            if (state.feasible[i].empty()) {
                return false;  // Wipeout
            }
        }
    }
    return true;
}

// Enforce target-weight uniqueness: when edge i with pair_sum P_i is assigned label
// L, the weight W = P_i + L is now claimed. For every other unassigned edge j with
// pair_sum P_j, the label (W - P_j) -- if it would also produce W -- must be
// removed from j's domain.
static bool forward_check_target_uniqueness(
    int assigned_edge,
    int assigned_label,
    int assigned_pair_sum,
    const std::vector<int>& pair_sums,
    const std::vector<bool>& edge_assigned,
    std::vector<std::set<int>>& feasible,
    std::vector<std::pair<int,int>>& removals) {

    int produced_weight = assigned_pair_sum + assigned_label;
    for (size_t j = 0; j < feasible.size(); ++j) {
        if ((int)j == assigned_edge) continue;
        if (edge_assigned[j]) continue;
        int forbidden_label = produced_weight - pair_sums[j];
        // forbidden_label, if in j's domain, would make j produce the same weight.
        auto it = feasible[j].find(forbidden_label);
        if (it != feasible[j].end()) {
            feasible[j].erase(it);
            removals.push_back({(int)j, forbidden_label});
            if (feasible[j].empty()) return false;
        }
    }
    return true;
}

void Solver::undo_forward_check(ConstraintState& state,
                                const std::vector<std::pair<int,int>>& removals) const {
    for (const auto& r : removals) {
        state.feasible[r.first].insert(r.second);
    }
}

// ============================================================
// LAYER 3: Greedy Initialization
// ============================================================

std::optional<std::vector<int>> Solver::greedy_assign(
    const std::vector<int>& pair_sums,
    const ConstraintState& initial_state) const {

    std::vector<int> assignment(pair_sums.size(), -1);
    std::unordered_set<int> used_labels;
    std::unordered_set<int> used_weights;

    // Walk edges in order; for each, pick smallest feasible label whose label and
    // resulting weight are both still available.
    for (size_t i = 0; i < pair_sums.size(); ++i) {
        int chosen = -1;
        for (int label : initial_state.feasible[i]) {
            if (used_labels.count(label)) continue;
            int weight = pair_sums[i] + label;
            if (used_weights.count(weight)) continue;
            chosen = label;
            break;
        }
        if (chosen < 0) return std::nullopt;  // Stuck

        assignment[i] = chosen;
        used_labels.insert(chosen);
        used_weights.insert(pair_sums[i] + chosen);
    }

    // Verify weights form an AP (they should automatically if both label and
    // weight uniqueness held and feasibility was correct, but double-check).
    std::vector<int> weights(pair_sums.size());
    for (size_t i = 0; i < pair_sums.size(); ++i) {
        weights[i] = pair_sums[i] + assignment[i];
    }
    std::sort(weights.begin(), weights.end());
    for (size_t i = 1; i < weights.size(); ++i) {
        if (weights[i] != weights[i-1] + 1) return std::nullopt;
    }

    return assignment;
}

// ============================================================
// LAYERS 4 & 5: Iterative Deepening DFS + MCV
// ============================================================

int Solver::select_next_edge_mcv(const ConstraintState& state,
                                 const std::vector<bool>& edge_assigned) const {
    int best = -1;
    size_t best_size = SIZE_MAX;
    for (size_t i = 0; i < state.feasible.size(); ++i) {
        if (edge_assigned[i]) continue;
        size_t s = state.feasible[i].size();
        if (s < best_size) {
            best_size = s;
            best = (int)i;
            if (best_size == 1) break;  // Can't do better
        }
    }
    return best;
}

Solver::DfsOutcome Solver::dfs(int depth,
                                int depth_limit,
                                const std::vector<int>& pair_sums,
                                ConstraintState& state,
                                std::vector<int>& assignment,
                                std::vector<bool>& edge_assigned,
                                const std::vector<int>* fixed_order) {

    // Check time budget
    if (time_exceeded()) return DfsOutcome::Failure;

    // Goal: every edge assigned
    int total = (int)pair_sums.size();
    int done = std::count(edge_assigned.begin(), edge_assigned.end(), true);
    if (done == total) return DfsOutcome::Success;

    // Iterative deepening cap
    if (depth >= depth_limit) return DfsOutcome::DepthExceeded;

    // Layer 6: memo check
    if (cfg_.use_memoization && memo_check_failure(done, assignment)) {
        memo_hits_count_++;
        return DfsOutcome::Failure;
    }

    // Pick next edge: either MCV or fixed order (for diversification)
    int edge_idx;
    if (fixed_order && (size_t)done < fixed_order->size()) {
        // Find next unassigned edge in the given order
        edge_idx = -1;
        for (int candidate : *fixed_order) {
            if (!edge_assigned[candidate]) {
                edge_idx = candidate;
                break;
            }
        }
        if (edge_idx < 0) return DfsOutcome::Failure;
    } else if (cfg_.use_mcv_ordering) {
        edge_idx = select_next_edge_mcv(state, edge_assigned);
        if (edge_idx < 0) return DfsOutcome::Failure;
    } else {
        // Default: leftmost unassigned
        edge_idx = -1;
        for (size_t i = 0; i < edge_assigned.size(); ++i) {
            if (!edge_assigned[i]) { edge_idx = (int)i; break; }
        }
        if (edge_idx < 0) return DfsOutcome::Failure;
    }

    // Snapshot the feasible set for this edge (we'll mutate it during forward checking)
    std::vector<int> options(state.feasible[edge_idx].begin(),
                             state.feasible[edge_idx].end());

    // Layer 8: Least-Constraining Value ordering.
    // For each candidate label, count how many other edges' domains it would shrink.
    // Try the LEAST constraining label first - it leaves the search most flexible,
    // increasing the chance that the first attempt succeeds.
    if (cfg_.use_lcv_value_ordering && options.size() > 1) {
        std::vector<std::pair<int,int>> scored;  // (constraint_count, label)
        scored.reserve(options.size());
        for (int label : options) {
            int weight = pair_sums[edge_idx] + label;
            int cost = 0;
            for (size_t j = 0; j < state.feasible.size(); ++j) {
                if ((int)j == edge_idx) continue;
                if (edge_assigned[j]) continue;
                // Would removing `label` shrink j's domain?
                if (state.feasible[j].count(label)) cost++;
                // Would removing the weight-clashing label shrink j's domain?
                int clash = weight - pair_sums[j];
                if (clash != label && state.feasible[j].count(clash)) cost++;
            }
            scored.push_back({cost, label});
        }
        std::sort(scored.begin(), scored.end());  // ascending by cost = least-constraining first
        options.clear();
        for (auto& p : scored) options.push_back(p.second);
    }

    bool any_depth_exceeded = false;

    for (int label : options) {
        // Tentatively assign
        assignment[edge_idx] = label;
        edge_assigned[edge_idx] = true;

        // Save the original domain of this edge so we can restore it
        std::set<int> saved_domain = state.feasible[edge_idx];
        state.feasible[edge_idx].clear();
        state.feasible[edge_idx].insert(label);

        // Forward check (a): remove `label` from other edges' domains
        std::vector<std::pair<int,int>> removals;
        bool fc_ok = forward_check(label, edge_assigned, state, removals);

        // Forward check (b): remove labels that would produce the same target weight
        if (fc_ok) {
            fc_ok = forward_check_target_uniqueness(
                edge_idx, label, pair_sums[edge_idx],
                pair_sums, edge_assigned, state.feasible, removals);
        }

        if (fc_ok) {
            DfsOutcome res = dfs(depth + 1, depth_limit, pair_sums,
                                 state, assignment, edge_assigned, fixed_order);
            if (res == DfsOutcome::Success) return res;
            if (res == DfsOutcome::DepthExceeded) any_depth_exceeded = true;
        }

        // Undo
        undo_forward_check(state, removals);
        state.feasible[edge_idx] = saved_domain;
        edge_assigned[edge_idx] = false;
        assignment[edge_idx] = -1;

        if (time_exceeded()) return DfsOutcome::Failure;
    }

    // No option worked. Record memo failure (only if we fully exhausted, not depth-cut).
    if (!any_depth_exceeded && cfg_.use_memoization) {
        memo_record_failure(done, assignment);
    }

    return any_depth_exceeded ? DfsOutcome::DepthExceeded : DfsOutcome::Failure;
}

bool Solver::iterative_deepening_search(
    const std::vector<int>& pair_sums,
    ConstraintState& state,
    std::vector<int>& assignment,
    const std::vector<int>& edge_order_hint) {

    int total = (int)pair_sums.size();
    std::vector<bool> edge_assigned(total, false);
    assignment.assign(total, -1);

    if (!cfg_.use_iterative_deepening) {
        // Single shot, depth = total
        DfsOutcome res = dfs(0, total, pair_sums, state, assignment, edge_assigned,
                             edge_order_hint.empty() ? nullptr : &edge_order_hint);
        return res == DfsOutcome::Success;
    }

    // For SEATL specifically: a complete solution requires exactly `total` edge
    // assignments (depth = total). Shallow attempts can never find a solution; they
    // would only waste time enumerating partial assignments. So instead of true
    // iterative deepening over solution depth, we use the depth limit as a *safety*
    // cap that protects against pathological recursion. The "iteration" only matters
    // if hard_cap < total, in which case no solution is reachable anyway.
    int hard_cap = std::min(cfg_.max_dfs_depth, total);

    std::fill(edge_assigned.begin(), edge_assigned.end(), false);
    std::fill(assignment.begin(), assignment.end(), -1);
    memo_failures_.clear();

    DfsOutcome res = dfs(0, hard_cap, pair_sums, state, assignment, edge_assigned,
                         edge_order_hint.empty() ? nullptr : &edge_order_hint);
    return res == DfsOutcome::Success;
}

// ============================================================
// DFS enumeration variant: records every solution found.
// Returns true iff search completed without being capped or timed out.
// ============================================================
bool Solver::dfs_enumerate(int depth,
                           int /*depth_limit*/,
                           const std::vector<int>& pair_sums,
                           const std::vector<int>& M,
                           const std::vector<int>& N,
                           int a,
                           ConstraintState& state,
                           std::vector<int>& assignment,
                           std::vector<bool>& edge_assigned,
                           std::vector<Solution>& solutions_out,
                           int max_solutions) {

    if (time_exceeded()) return false;
    if ((int)solutions_out.size() >= max_solutions) return false;

    int total = (int)pair_sums.size();
    int done  = std::count(edge_assigned.begin(), edge_assigned.end(), true);

    // Goal reached: record the solution and continue searching for more.
    if (done == total) {
        Solution sol;
        sol.a = a;
        sol.set_m = M;
        sol.set_n = N;
        int idx = 0;
        for (int u : M) {
            for (int v : N) {
                sol.edge_labeling.push_back({{u, v}, assignment[idx]});
                sol.edge_weights.push_back(u + v + assignment[idx]);
                idx++;
            }
        }
        std::sort(sol.edge_weights.begin(), sol.edge_weights.end());
        solutions_out.push_back(std::move(sol));
        return true;  // Continue (caller will keep exploring siblings)
    }

    // Pick next edge using MCV when enabled.
    int edge_idx = -1;
    if (cfg_.use_mcv_ordering) {
        edge_idx = select_next_edge_mcv(state, edge_assigned);
    } else {
        for (int i = 0; i < total; ++i) {
            if (!edge_assigned[i]) { edge_idx = i; break; }
        }
    }
    if (edge_idx < 0) return true;

    std::vector<int> options(state.feasible[edge_idx].begin(),
                             state.feasible[edge_idx].end());

    for (int label : options) {
        if (time_exceeded()) return false;
        if ((int)solutions_out.size() >= max_solutions) return false;

        assignment[edge_idx] = label;
        edge_assigned[edge_idx] = true;

        std::set<int> saved_domain = state.feasible[edge_idx];
        state.feasible[edge_idx].clear();
        state.feasible[edge_idx].insert(label);

        std::vector<std::pair<int,int>> removals;
        bool fc_ok = forward_check(label, edge_assigned, state, removals);
        if (fc_ok) {
            fc_ok = forward_check_target_uniqueness(
                edge_idx, label, pair_sums[edge_idx],
                pair_sums, edge_assigned, state.feasible, removals);
        }

        if (fc_ok) {
            dfs_enumerate(depth + 1, total, pair_sums, M, N, a,
                          state, assignment, edge_assigned,
                          solutions_out, max_solutions);
        }

        // Always undo and try the next label - we want all solutions.
        undo_forward_check(state, removals);
        state.feasible[edge_idx] = saved_domain;
        edge_assigned[edge_idx] = false;
        assignment[edge_idx] = -1;
    }
    return true;
}

// ============================================================
// LAYER 6: Memoization
// ============================================================

unsigned long long Solver::memo_key(int edge_idx,
                                    const std::vector<int>& assignment) const {
    // Hash combines: edge_idx + the multiset of used labels.
    // Since each label is unique, we can XOR-hash them with mixing.
    unsigned long long h = (unsigned long long)edge_idx * 0x9E3779B97F4A7C15ULL;
    for (int v : assignment) {
        if (v < 0) continue;
        unsigned long long x = (unsigned long long)v;
        // Splitmix-style mixing
        x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
        x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
        x = x ^ (x >> 31);
        h ^= x;
    }
    return h;
}

bool Solver::memo_check_failure(int edge_idx,
                                const std::vector<int>& assignment) const {
    unsigned long long k = memo_key(edge_idx, assignment);
    return memo_failures_.count(k) > 0;
}

void Solver::memo_record_failure(int edge_idx,
                                 const std::vector<int>& assignment) {
    if ((int)memo_failures_.size() >= cfg_.memoization_limit) return;
    memo_failures_.insert(memo_key(edge_idx, assignment));
}

// ============================================================
// LAYER 7: Restart with diversification
// ============================================================

bool Solver::solve_with_restarts(const std::vector<int>& pair_sums,
                                 ConstraintState& state,
                                 std::vector<int>& assignment,
                                 int& restarts_used) {

    int total = (int)pair_sums.size();

    // Attempt 1: standard MCV ordering, no fixed order.
    {
        ConstraintState working = state;
        std::vector<int> empty_order;
        if (iterative_deepening_search(pair_sums, working, assignment, empty_order)) {
            return true;
        }
    }

    if (!cfg_.use_restart) return false;

    // Subsequent attempts: shuffled fixed orders.
    std::mt19937 rng(0xC0FFEE);
    for (int attempt = 1; attempt < cfg_.restart_attempts; ++attempt) {
        if (time_exceeded()) return false;
        restarts_used++;

        std::vector<int> order(total);
        std::iota(order.begin(), order.end(), 0);
        std::shuffle(order.begin(), order.end(), rng);

        ConstraintState working = state;
        if (iterative_deepening_search(pair_sums, working, assignment, order)) {
            return true;
        }
    }
    return false;
}

// ============================================================
// Helpers: pair sums, targets, partitions, result building
// ============================================================

std::vector<int> Solver::compute_pair_sums(const std::vector<int>& M,
                                           const std::vector<int>& N) const {
    std::vector<int> ps;
    ps.reserve(mn_);
    for (int u : M) {
        for (int v : N) {
            ps.push_back(u + v);
        }
    }
    return ps;
}

std::vector<int> Solver::compute_targets(int a) const {
    std::vector<int> t(mn_);
    for (int i = 0; i < mn_; ++i) t[i] = a + i;
    return t;
}

void Solver::generate_partitions(
    std::vector<std::pair<std::vector<int>, std::vector<int>>>& out) const {

    int total = m_ + n_;
    std::vector<int> all(total);
    std::iota(all.begin(), all.end(), 1);

    std::vector<int> mask(total, 0);
    for (int i = 0; i < m_; ++i) mask[i] = 1;
    std::sort(mask.begin(), mask.end());

    do {
        std::vector<int> M, N;
        for (int i = 0; i < total; ++i) {
            if (mask[i]) M.push_back(all[i]);
            else N.push_back(all[i]);
        }
        // Canonical: only keep partitions where M contains label 1.
        if (m_ == n_) {
            if (std::find(M.begin(), M.end(), 1) == M.end()) continue;
        }
        out.push_back({M, N});
    } while (std::next_permutation(mask.begin(), mask.end()));

    // Layer 10: order partitions by promise.
    // Partitions where pair sums span a wide range tend to be solvable because
    // they can hit a wider variety of target weights. Score each partition by
    // the variance/spread of its pair sums.
    if (cfg_.use_partition_ordering) {
        struct Scored { int score; size_t index; };
        std::vector<Scored> scored;
        scored.reserve(out.size());
        for (size_t i = 0; i < out.size(); ++i) {
            const auto& M = out[i].first;
            const auto& N = out[i].second;
            // Compute spread of pair sums: max - min
            int min_ps = INT_MAX, max_ps = INT_MIN;
            for (int u : M) for (int v : N) {
                int ps = u + v;
                if (ps < min_ps) min_ps = ps;
                if (ps > max_ps) max_ps = ps;
            }
            // We want partitions whose pair-sum range matches mn (so each weight
            // is reachable by some edge). The "ideal" spread is mn-1. Penalize
            // partitions that deviate from that.
            int spread = max_ps - min_ps;
            int ideal = mn_ - 1;
            scored.push_back({std::abs(spread - ideal), i});
        }
        std::sort(scored.begin(), scored.end(),
                  [](const Scored& a, const Scored& b){ return a.score < b.score; });
        std::vector<std::pair<std::vector<int>, std::vector<int>>> reordered;
        reordered.reserve(out.size());
        for (auto& s : scored) reordered.push_back(std::move(out[s.index]));
        out = std::move(reordered);
    }
}

Result Solver::build_result(const std::vector<int>& M,
                            const std::vector<int>& N,
                            const std::vector<int>& assignment,
                            int a) const {
    Result r;
    r.found = true;
    r.a = a;
    r.set_m = M;
    r.set_n = N;

    // Reconstruct edges in the same order compute_pair_sums uses
    int idx = 0;
    for (int u : M) {
        for (int v : N) {
            r.edge_labeling.push_back({{u, v}, assignment[idx]});
            r.edge_weights.push_back(u + v + assignment[idx]);
            idx++;
        }
    }
    std::sort(r.edge_weights.begin(), r.edge_weights.end());
    return r;
}

// ============================================================
// MAIN ENTRY: solve()
// ============================================================

Result Solver::solve() {
    Result result;
    start_time_ = std::chrono::steady_clock::now();
    memo_failures_.clear();
    memo_hits_count_ = 0;

    // Generate all partitions
    std::vector<std::pair<std::vector<int>, std::vector<int>>> partitions;
    generate_partitions(partitions);

    if (cfg_.verbose) {
        std::cout << "[SEATL] K_{" << m_ << "," << n_ << "}: "
                  << partitions.size() << " canonical partitions to test\n";
    }

    int restarts_used = 0;

    for (auto& [M, N] : partitions) {
        if (total_time_exceeded()) {
            result.timed_out = true;
            result.status = "Total time limit exceeded";
            break;
        }

        // Set a per-partition deadline so one bad partition can't burn the whole budget.
        has_partition_deadline_ = true;
        partition_deadline_ = std::chrono::steady_clock::now() +
            std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                std::chrono::duration<double>(cfg_.per_partition_time));

        result.partitions_tested++;

        // ----- LAYER 1: Master formula -----
        int a = -1;
        if (cfg_.use_master_formula) {
            a = compute_starting_weight(M, N);
            if (a < 0) {
                result.partitions_pruned_layer1++;
                continue;
            }
        } else {
            // Without master formula, we still need 'a' to define targets.
            // Compute it anyway (it's free); just don't reject on non-integer.
            a = compute_starting_weight(M, N);
            if (a < 0) continue;  // Non-integer -> truly impossible regardless
        }

        std::vector<int> pair_sums = compute_pair_sums(M, N);
        std::vector<int> targets   = compute_targets(a);

        // ----- LAYER 2: AC-3 -----
        ConstraintState state;
        if (cfg_.use_ac3) {
            if (!build_initial_feasibility(pair_sums, targets, state)) {
                result.partitions_pruned_layer2++;
                continue;
            }
            if (!ac3_propagate(pair_sums, targets, state)) {
                result.partitions_pruned_layer2++;
                continue;
            }
        } else {
            // Build the table without propagation.
            if (!build_initial_feasibility(pair_sums, targets, state)) {
                result.partitions_pruned_layer2++;
                continue;
            }
        }

        // ----- LAYER 3: Greedy -----
        if (cfg_.use_greedy) {
            auto greedy = greedy_assign(pair_sums, state);
            if (greedy) {
                Result solution = build_result(M, N, *greedy, a);
                solution.partitions_tested        = result.partitions_tested;
                solution.partitions_pruned_layer1 = result.partitions_pruned_layer1;
                solution.partitions_pruned_layer2 = result.partitions_pruned_layer2;
                solution.greedy_successes         = result.greedy_successes + 1;
                solution.dfs_invocations          = result.dfs_invocations;
                solution.status = "Solved by greedy";
                result = solution;
                goto done;
            }
        }

        // ----- LAYER 11: Local search BEFORE expensive DFS -----
        // For square graphs at scale, DFS often grinds. Local search can find
        // a solution in hundreds of milliseconds when a solution exists.
        if (cfg_.use_local_search) {
            std::vector<int> assignment;
            if (local_search(pair_sums, targets, state, assignment)) {
                Result solution = build_result(M, N, assignment, a);
                solution.partitions_tested        = result.partitions_tested;
                solution.partitions_pruned_layer1 = result.partitions_pruned_layer1;
                solution.partitions_pruned_layer2 = result.partitions_pruned_layer2;
                solution.greedy_successes         = result.greedy_successes;
                solution.dfs_invocations          = result.dfs_invocations;
                solution.status = "Solved by local search";
                result = solution;
                goto done;
            }
        }

        // ----- LAYERS 4-7: Search -----
        {
            std::vector<int> assignment;
            result.dfs_invocations++;
            bool ok = solve_with_restarts(pair_sums, state, assignment, restarts_used);
            if (ok) {
                Result solution = build_result(M, N, assignment, a);
                solution.partitions_tested        = result.partitions_tested;
                solution.partitions_pruned_layer1 = result.partitions_pruned_layer1;
                solution.partitions_pruned_layer2 = result.partitions_pruned_layer2;
                solution.greedy_successes         = result.greedy_successes;
                solution.dfs_invocations          = result.dfs_invocations;
                solution.status = "Solved by DFS";
                result = solution;
                goto done;
            }
        }
    }

    if (!result.found && !result.timed_out) {
        result.status = "No solution found across all partitions";
    }

done:
    has_partition_deadline_ = false;
    auto end = std::chrono::steady_clock::now();
    result.elapsed_ms = std::chrono::duration<double, std::milli>(end - start_time_).count();
    result.memo_hits = memo_hits_count_;
    result.restarts_used = restarts_used;
    return result;
}

// ============================================================
// solve_all() - enumerate as many SEATL solutions as we can
// ============================================================
//
// Strategy: run a systematic DFS-only enumeration over every canonical partition.
// Local search is intentionally NOT used here - it samples solutions but cannot
// guarantee enumeration. The trade-off: enumeration is much slower than finding
// just one solution, especially for large graphs. Use max_solutions and
// time_limit_seconds to bound the work.
//
// Completeness: the result is "complete" only if (a) no time limit was hit, AND
// (b) max_solutions cap was not reached, AND (c) every partition was processed.

EnumerationResult Solver::solve_all() {
    EnumerationResult result;
    start_time_ = std::chrono::steady_clock::now();
    memo_failures_.clear();
    memo_hits_count_ = 0;

    std::vector<std::pair<std::vector<int>, std::vector<int>>> partitions;
    generate_partitions(partitions);

    if (cfg_.verbose) {
        std::cout << "[SEATL] enumerating K_{" << m_ << "," << n_ << "}: "
                  << partitions.size() << " canonical partitions, "
                  << "max_solutions=" << cfg_.max_solutions
                  << ", time_limit=" << cfg_.time_limit_seconds << "s\n";
    }

    bool fully_complete = true;

    for (auto& [M, N] : partitions) {
        if (total_time_exceeded()) {
            result.timed_out = true;
            fully_complete = false;
            break;
        }
        if ((int)result.solutions.size() >= cfg_.max_solutions) {
            result.capped = true;
            fully_complete = false;
            break;
        }

        // Per-partition deadline so one stubborn partition can't hog the budget.
        has_partition_deadline_ = true;
        partition_deadline_ = std::chrono::steady_clock::now() +
            std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                std::chrono::duration<double>(cfg_.per_partition_time));

        result.partitions_tested++;

        // Layer 1: master formula
        int a = compute_starting_weight(M, N);
        if (a < 0) {
            result.partitions_pruned++;
            continue;
        }

        std::vector<int> pair_sums = compute_pair_sums(M, N);
        std::vector<int> targets   = compute_targets(a);

        // Layer 2: AC-3
        ConstraintState state;
        if (!build_initial_feasibility(pair_sums, targets, state)) {
            result.partitions_pruned++;
            continue;
        }
        if (cfg_.use_ac3 && !ac3_propagate(pair_sums, targets, state)) {
            result.partitions_pruned++;
            continue;
        }

        // Enumerate via DFS
        int total = (int)pair_sums.size();
        std::vector<int> assignment(total, -1);
        std::vector<bool> edge_assigned(total, false);

        size_t before = result.solutions.size();
        bool exhaustive = dfs_enumerate(0, total, pair_sums, M, N, a,
                                         state, assignment, edge_assigned,
                                         result.solutions, cfg_.max_solutions);
        if (!exhaustive) {
            // Either time exceeded or cap hit. We've recorded what we found.
            fully_complete = false;
        }
        if ((int)result.solutions.size() >= cfg_.max_solutions) {
            result.capped = true;
            fully_complete = false;
            break;
        }

        if (cfg_.verbose) {
            size_t added = result.solutions.size() - before;
            if (added > 0) {
                std::cout << "  partition " << result.partitions_tested
                          << " (M=" << M.size() << ",N=" << N.size() << "): "
                          << added << " solutions found\n";
            }
        }
    }

    has_partition_deadline_ = false;
    result.complete = fully_complete && !result.timed_out && !result.capped;

    auto end = std::chrono::steady_clock::now();
    result.elapsed_ms = std::chrono::duration<double, std::milli>(end - start_time_).count();

    if (result.complete) {
        result.status = "Complete enumeration: " + std::to_string(result.solutions.size())
                      + " solution(s) found";
    } else if (result.capped) {
        result.status = "Stopped at max_solutions cap (" + std::to_string(cfg_.max_solutions) + ")";
    } else if (result.timed_out) {
        result.status = "Time limit exceeded; partial enumeration";
    } else {
        result.status = "Partial enumeration (some partitions had per-partition timeout)";
    }
    return result;
}

// ============================================================
// LAYER 11: Local search (min-conflicts repair)
// ============================================================
//
// Reframe SEATL as: find a permutation pi of {0..mn-1} such that for each edge i,
// edge i is assigned label = (a + pi(i)) - pair_sums[i], and all those labels are
// distinct integers in the valid edge-label range. The labels being "in range" is
// already guaranteed by the feasibility table; the labels being "distinct" is the
// hard global constraint.
//
// Min-conflicts:
//   1. Start with a random permutation.
//   2. Compute conflicts: edges whose computed labels collide with another edge's.
//   3. Pick a conflicted edge; swap its target with another edge such that the
//      number of conflicts strictly decreases (or, with small probability, accept
//      a non-improving swap to escape local minima).
//   4. Repeat until either no conflicts remain (solution!) or iteration cap hit.

bool Solver::local_search(const std::vector<int>& pair_sums,
                          const std::vector<int>& targets,
                          const ConstraintState& state,
                          std::vector<int>& assignment_out) {

    int n = (int)pair_sums.size();
    std::mt19937 rng(0xDEADBEEF);

    // Pre-compute: for each edge, which target indices yield a label in its feasible set?
    std::vector<std::vector<int>> allowed_targets(n);  // edge -> list of target indices
    for (int i = 0; i < n; ++i) {
        for (int t = 0; t < n; ++t) {
            int needed = targets[t] - pair_sums[i];
            if (state.feasible[i].count(needed)) {
                allowed_targets[i].push_back(t);
            }
        }
        if (allowed_targets[i].empty()) return false;  // Infeasible
    }

    auto compute_label = [&](int edge_i, int target_idx) {
        return targets[target_idx] - pair_sums[edge_i];
    };

    auto count_label_conflicts = [&](const std::vector<int>& perm) {
        // perm[i] = target index assigned to edge i. Count edges whose label collides.
        std::unordered_map<int,int> label_count;
        for (int i = 0; i < n; ++i) {
            int lbl = compute_label(i, perm[i]);
            label_count[lbl]++;
        }
        int conflicts = 0;
        for (int i = 0; i < n; ++i) {
            int lbl = compute_label(i, perm[i]);
            if (label_count[lbl] > 1) conflicts++;
        }
        return conflicts;
    };

    for (int attempt = 0; attempt < cfg_.local_search_attempts; ++attempt) {
        if (time_exceeded()) return false;

        // Initialize with a random valid permutation.
        // Build a random feasible assignment using a try-and-retry approach.
        std::vector<int> perm(n, -1);
        std::vector<bool> target_used(n, false);
        bool init_ok = true;
        std::vector<int> order(n);
        std::iota(order.begin(), order.end(), 0);
        std::shuffle(order.begin(), order.end(), rng);
        for (int i : order) {
            std::vector<int> shuffled_targets = allowed_targets[i];
            std::shuffle(shuffled_targets.begin(), shuffled_targets.end(), rng);
            int chosen = -1;
            for (int t : shuffled_targets) {
                if (!target_used[t]) { chosen = t; break; }
            }
            if (chosen < 0) { init_ok = false; break; }
            perm[i] = chosen;
            target_used[chosen] = true;
        }
        if (!init_ok) continue;

        // Min-conflicts iterative repair
        for (int iter = 0; iter < cfg_.local_search_iterations; ++iter) {
            if (iter % 100 == 0 && time_exceeded()) return false;

            // Compute conflicts: which edges share a label with another edge?
            std::unordered_map<int,std::vector<int>> label_to_edges;
            for (int i = 0; i < n; ++i) {
                int lbl = compute_label(i, perm[i]);
                label_to_edges[lbl].push_back(i);
            }
            std::vector<int> conflicted;
            for (auto& kv : label_to_edges) {
                if (kv.second.size() > 1) {
                    for (int e : kv.second) conflicted.push_back(e);
                }
            }
            if (conflicted.empty()) {
                // SOLUTION FOUND. Convert perm to label assignment.
                assignment_out.assign(n, 0);
                for (int i = 0; i < n; ++i) {
                    assignment_out[i] = compute_label(i, perm[i]);
                }
                return true;
            }

            // Pick a random conflicted edge
            int e1 = conflicted[std::uniform_int_distribution<int>(0, (int)conflicted.size()-1)(rng)];

            // Try to find a swap partner e2 such that swapping perm[e1] and perm[e2]
            // reduces total conflicts. Search over all edges; accept first improvement.
            int current_conflicts = (int)conflicted.size();
            int best_delta = 0;
            int best_partner = -1;
            for (int e2 = 0; e2 < n; ++e2) {
                if (e2 == e1) continue;
                // Check that the swap is feasible: e2 can take perm[e1]'s target,
                // and e1 can take perm[e2]'s target.
                int t1 = perm[e1], t2 = perm[e2];
                bool e1_can_take_t2 = false, e2_can_take_t1 = false;
                for (int t : allowed_targets[e1]) if (t == t2) { e1_can_take_t2 = true; break; }
                if (!e1_can_take_t2) continue;
                for (int t : allowed_targets[e2]) if (t == t1) { e2_can_take_t1 = true; break; }
                if (!e2_can_take_t1) continue;

                // Try the swap, count new conflicts, undo if no improvement.
                std::swap(perm[e1], perm[e2]);
                int new_conflicts = count_label_conflicts(perm);
                int delta = current_conflicts - new_conflicts;
                std::swap(perm[e1], perm[e2]);  // undo
                if (delta > best_delta) {
                    best_delta = delta;
                    best_partner = e2;
                }
            }

            if (best_partner >= 0) {
                std::swap(perm[e1], perm[best_partner]);
            } else {
                // Local minimum - accept a random swap to escape (5% chance)
                std::uniform_int_distribution<int> coin(0, 19);
                if (coin(rng) == 0) {
                    int e2 = std::uniform_int_distribution<int>(0, n-1)(rng);
                    if (e2 != e1) {
                        int t1 = perm[e1], t2 = perm[e2];
                        bool e1_ok = false, e2_ok = false;
                        for (int t : allowed_targets[e1]) if (t == t2) { e1_ok = true; break; }
                        for (int t : allowed_targets[e2]) if (t == t1) { e2_ok = true; break; }
                        if (e1_ok && e2_ok) std::swap(perm[e1], perm[e2]);
                    }
                }
            }
        }
    }
    return false;
}

} // namespace seatl
