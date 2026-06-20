/*
 * seatl_solver_xl.cpp - Implementation of the XL (Layer 12) solver.
 *
 * Key optimizations for K_{20,20}+:
 *  1. Skip partition enumeration (which has 10^11 candidates)
 *  2. Use bitset-based fast lookups instead of std::set
 *  3. Parallel tabu search via OpenMP
 *  4. Smart partition seeding (start with the "block" partition that has
 *     the densest pair-sum distribution)
 */

#include "seatl_solver_xl.hpp"

#include <algorithm>
#include <numeric>
#include <random>
#include <chrono>
#include <iostream>
#include <vector>
#include <cstring>
#include <climits>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace seatl {

// ============================================================
// Constructor
// ============================================================
SolverXL::SolverXL(int m, int n, const ConfigXL& cfg)
    : m_(m), n_(n), mn_(m * n),
      min_label_(m + n + 1), max_label_(m + n + m * n),
      cfg_(cfg) {}

// ============================================================
// Time check
// ============================================================
bool SolverXL::time_exceeded() const {
    auto now = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration<double>(now - start_time_).count();
    return elapsed >= cfg_.time_limit_seconds;
}

// ============================================================
// Master Formula (Layer 1)
// ============================================================
int SolverXL::compute_a(const std::vector<int>& M, const std::vector<int>& N) const {
    long long sum_m = std::accumulate(M.begin(), M.end(), 0LL);
    long long sum_n = std::accumulate(N.begin(), N.end(), 0LL);
    long long e_lo = min_label_, e_hi = max_label_;
    long long sum_e = (e_lo + e_hi) * mn_ / 2;
    long long num = (long long)n_ * sum_m
                  + (long long)m_ * sum_n
                  + sum_e
                  - (long long)mn_ * (mn_ - 1) / 2;
    if (num % mn_ != 0) return -1;
    long long a = num / mn_;
    if (a < 1) return -1;
    return (int)a;
}

// ============================================================
// Pair sums
// ============================================================
std::vector<int> SolverXL::compute_pair_sums(const std::vector<int>& M,
                                              const std::vector<int>& N) const {
    std::vector<int> ps;
    ps.reserve(mn_);
    for (int u : M) for (int v : N) ps.push_back(u + v);
    return ps;
}

// ============================================================
// Smart partition seeding
// ============================================================
// For K_{m,n} large, we don't enumerate all C(m+n, m) partitions.
// Instead, we try the "block" partition first, then jitter.
//
// For rectangular K_{m,n} (m != n), the simple seeds may not yield
// integer 'a'. We then search for valid seeds by jittering until we
// find some.
// ============================================================
void SolverXL::generate_seed_partitions(
    std::vector<std::pair<std::vector<int>, std::vector<int>>>& out) const {

    std::vector<int> all(m_ + n_);
    std::iota(all.begin(), all.end(), 1);

    auto try_add = [&](std::vector<int> M, std::vector<int> N) -> bool {
        std::sort(M.begin(), M.end());
        std::sort(N.begin(), N.end());
        // Dedupe
        for (auto& [eM, eN] : out) if (eM == M && eN == N) return false;
        if (compute_a(M, N) >= 1) {
            out.push_back({M, N});
            return true;
        }
        return false;
    };

    // SEED 1: Block partition M={1..m}, N={m+1..m+n}
    {
        std::vector<int> M(all.begin(), all.begin() + m_);
        std::vector<int> N(all.begin() + m_, all.end());
        try_add(M, N);
    }

    // SEED 2: Interleaved partition M={1,3,5,...}, N={2,4,6,...}
    {
        std::vector<int> M, N;
        for (int i = 1; i <= m_ + n_; ++i) {
            if (M.size() < (size_t)m_ && (i % 2 == 1 || N.size() == (size_t)n_)) {
                M.push_back(i);
            } else {
                N.push_back(i);
            }
        }
        try_add(M, N);
    }

    // SEED 3: Reverse block partition (only if m == n, otherwise sizes don't match)
    if (m_ == n_) {
        std::vector<int> M(all.begin() + n_, all.end());
        std::vector<int> N(all.begin(), all.begin() + n_);
        try_add(M, N);
    }

    // SEED 4: Random search for valid partitions
    // For rectangular graphs, simple seeds may not pass the master formula.
    // We jitter starting from the block partition until we find valid 'a'.
    std::mt19937 rng(0x9E3779B9);
    int max_attempts = std::max(500, m_ * n_);
    int target_count = std::max(8, cfg_.partition_jitter);

    // Start from block partition
    std::vector<int> M(all.begin(), all.begin() + m_);
    std::vector<int> N(all.begin() + m_, all.end());

    for (int attempt = 0; attempt < max_attempts && (int)out.size() < target_count; ++attempt) {
        // Swap one random element between M and N
        std::uniform_int_distribution<int> di(0, m_ - 1), dj(0, n_ - 1);
        int i = di(rng), j = dj(rng);
        std::swap(M[i], N[j]);
        try_add(M, N);
    }

    // If still nothing, try harder with multi-swap jitter
    if (out.empty()) {
        for (int attempt = 0; attempt < max_attempts * 2; ++attempt) {
            M.assign(all.begin(), all.begin() + m_);
            N.assign(all.begin() + m_, all.end());
            // Multiple random swaps
            std::uniform_int_distribution<int> di(0, m_ - 1), dj(0, n_ - 1);
            int num_swaps = 1 + (attempt % 5);
            for (int k = 0; k < num_swaps; ++k) {
                int i = di(rng), j = dj(rng);
                std::swap(M[i], N[j]);
            }
            if (try_add(M, N)) {
                if ((int)out.size() >= target_count) break;
            }
        }
    }
}

// ============================================================
// Tabu Search Worker (Layer 12 core) — OPTIMIZED VERSION
// ============================================================
// Key optimization: O(1) incremental scoring when evaluating swaps.
//
// For a swap of edges (i, j):
//   Old weight of i = pair_sums[i] + assignment[i]
//   Old weight of j = pair_sums[j] + assignment[j]
//   New weight of i = pair_sums[i] + assignment[j]
//   New weight of j = pair_sums[j] + assignment[i]
//
// We track `slot_count[t]` = how many edges hit target slot t.
// Score = sum over slots of |count[t] - 1| (target should be hit exactly once).
//
// A swap only changes 4 slots: old_w_i, old_w_j, new_w_i, new_w_j.
// We can compute the delta in O(1) without rescanning.
// ============================================================
bool SolverXL::tabu_worker(const std::vector<int>& pair_sums,
                            int a,
                            std::vector<int>& best_assignment,
                            int seed,
                            std::chrono::steady_clock::time_point deadline) {
    int mn = mn_;
    int min_label = min_label_;

    std::mt19937 rng(seed);

    // Initialize: shuffle labels onto edges
    std::vector<int> assignment(mn);
    {
        std::vector<int> labels(mn);
        std::iota(labels.begin(), labels.end(), min_label);
        std::shuffle(labels.begin(), labels.end(), rng);
        for (int i = 0; i < mn; ++i) assignment[i] = labels[i];
    }

    // The "slot" range covers the target window [a, a+mn-1]
    // But weights can be out-of-range, so we use a wider buffer
    // weight = u + v + label where u,v ∈ [1, m+n], label ∈ [m+n+1, m+n+mn]
    // min weight = 1+2+(m+n+1) = m+n+4
    // max weight = (m+n-1)+(m+n)+(m+n+mn) = 3(m+n)+mn-1
    int min_possible_weight = 1 + 2 + min_label;
    int max_possible_weight = (m_ + n_ - 1) + (m_ + n_) + (m_ + n_ + mn_);
    int slot_offset = min_possible_weight; // index 0 corresponds to this weight
    int slot_count_size = max_possible_weight - min_possible_weight + 2;

    std::vector<int> slot_count(slot_count_size, 0);

    auto slot_idx = [&](int weight) -> int {
        int s = weight - slot_offset;
        if (s < 0) return 0;
        if (s >= slot_count_size) return slot_count_size - 1;
        return s;
    };

    // Target window in slot-space: [a-slot_offset, a-slot_offset+mn-1]
    int target_lo = a - slot_offset;
    int target_hi = a - slot_offset + mn - 1;

    // Initialize slot_count
    for (int i = 0; i < mn; ++i) {
        int w = pair_sums[i] + assignment[i];
        slot_count[slot_idx(w)]++;
    }

    // Compute current score from slot_count
    auto compute_score = [&]() -> int {
        int s = 0;
        // Target slots: want exactly 1
        for (int t = target_lo; t <= target_hi; ++t) {
            s += std::abs(slot_count[t] - 1);
        }
        // Non-target slots: want exactly 0
        for (int t = 0; t < target_lo; ++t) s += slot_count[t];
        for (int t = target_hi + 1; t < slot_count_size; ++t) s += slot_count[t];
        return s;
    };

    int current_score = compute_score();
    if (current_score == 0) {
        best_assignment = assignment;
        return true;
    }

    // Score-delta from updating one slot's count (before -> after)
    auto slot_contribution = [&](int t, int count) -> int {
        if (t >= target_lo && t <= target_hi) return std::abs(count - 1);
        else return count;
    };

    // O(1) delta computation for swapping edges i and j
    auto swap_delta = [&](int i, int j) -> int {
        if (i == j) return 0;
        int li = assignment[i], lj = assignment[j];
        int old_wi = pair_sums[i] + li;
        int old_wj = pair_sums[j] + lj;
        int new_wi = pair_sums[i] + lj;
        int new_wj = pair_sums[j] + li;

        int si_old = slot_idx(old_wi);
        int sj_old = slot_idx(old_wj);
        int si_new = slot_idx(new_wi);
        int sj_new = slot_idx(new_wj);

        // We're removing 1 from si_old, 1 from sj_old, then adding 1 to si_new and sj_new
        // But some slots might overlap — handle carefully

        // Collect unique slots and their before/after counts
        struct SlotChange { int slot; int delta; };
        SlotChange changes[4] = {
            {si_old, -1}, {sj_old, -1}, {si_new, +1}, {sj_new, +1}
        };

        // Compute delta in score
        int delta = 0;
        // For each affected slot, compute the change in its contribution
        // We need to handle the case where the same slot appears multiple times
        // Strategy: for each unique slot, compute total net change

        // Simple approach: iterate, tracking already-processed slots
        int processed[4] = {-1, -1, -1, -1};
        int np = 0;
        for (int c = 0; c < 4; ++c) {
            int sl = changes[c].slot;
            bool seen = false;
            for (int p = 0; p < np; ++p) {
                if (processed[p] == sl) { seen = true; break; }
            }
            if (seen) continue;
            processed[np++] = sl;

            // Compute net delta for this slot
            int net = 0;
            for (int c2 = 0; c2 < 4; ++c2) {
                if (changes[c2].slot == sl) net += changes[c2].delta;
            }
            if (net == 0) continue;

            int old_cnt = slot_count[sl];
            int new_cnt = old_cnt + net;
            delta += slot_contribution(sl, new_cnt) - slot_contribution(sl, old_cnt);
        }

        return delta;
    };

    // Apply swap: updates slot_count and assignment
    auto apply_swap = [&](int i, int j) {
        int li = assignment[i], lj = assignment[j];
        int old_wi = pair_sums[i] + li;
        int old_wj = pair_sums[j] + lj;
        int new_wi = pair_sums[i] + lj;
        int new_wj = pair_sums[j] + li;

        slot_count[slot_idx(old_wi)]--;
        slot_count[slot_idx(old_wj)]--;
        slot_count[slot_idx(new_wi)]++;
        slot_count[slot_idx(new_wj)]++;

        std::swap(assignment[i], assignment[j]);
    };

    std::vector<int> best_known = assignment;
    int best_score = current_score;

    // Tabu table: tabu_until[i*mn + j]
    std::vector<int> tabu_until(mn * mn, 0);
    auto tabu_set = [&](int i, int j, int until) {
        tabu_until[i * mn + j] = until;
        tabu_until[j * mn + i] = until;
    };
    auto tabu_get = [&](int i, int j) -> int {
        return tabu_until[i * mn + j];
    };

    int iter = 0;
    int stagnation = 0;
    const int max_stagnation = std::max(2000, mn * 10);
    const int sample_size = std::max(100, mn);

    std::uniform_int_distribution<int> di(0, mn - 1);

    for (; iter < cfg_.tabu_iterations; ++iter) {
        // Check deadline periodically
        if ((iter & 4095) == 0) {
            if (std::chrono::steady_clock::now() >= deadline) break;
        }

        // Sample candidate swaps
        int best_i = -1, best_j = -1;
        int best_delta = INT_MAX;

        for (int s = 0; s < sample_size; ++s) {
            int i = di(rng), j = di(rng);
            if (i == j) continue;
            bool is_tabu = (tabu_get(i, j) > iter);

            int delta = swap_delta(i, j);

            // Aspiration: accept tabu move if it improves best known
            if (is_tabu && current_score + delta >= best_score) continue;

            if (delta < best_delta) {
                best_delta = delta;
                best_i = i;
                best_j = j;
                if (delta < -2) break;  // good enough
            }
        }

        if (best_i < 0) {
            stagnation++;
            if (stagnation > 50) {
                // Diversify: random swaps
                for (int s = 0; s < 3; ++s) {
                    int i = di(rng), j = di(rng);
                    if (i != j) apply_swap(i, j);
                }
                current_score = compute_score();
                stagnation = 0;
            }
            continue;
        }

        // Apply best swap
        apply_swap(best_i, best_j);
        current_score += best_delta;
        tabu_set(best_i, best_j, iter + cfg_.tabu_tenure);

        if (current_score < best_score) {
            best_score = current_score;
            best_known = assignment;
            stagnation = 0;
            if (best_score == 0) {
                best_assignment = best_known;
                return true;
            }
        } else {
            stagnation++;
            if (stagnation > max_stagnation) {
                // Restart with shuffle
                std::vector<int> labels(mn);
                std::iota(labels.begin(), labels.end(), min_label);
                std::shuffle(labels.begin(), labels.end(), rng);
                assignment = labels;
                std::fill(slot_count.begin(), slot_count.end(), 0);
                for (int i = 0; i < mn; ++i) {
                    slot_count[slot_idx(pair_sums[i] + assignment[i])]++;
                }
                current_score = compute_score();
                stagnation = 0;
                std::fill(tabu_until.begin(), tabu_until.end(), 0);
            }
        }
    }

    best_assignment = best_known;
    return best_score == 0;
}

// ============================================================
// Parallel tabu search
// ============================================================
bool SolverXL::tabu_search(const std::vector<int>& pair_sums,
                            int a,
                            std::vector<int>& assignment_out,
                            int thread_seed) {

    int num_threads = cfg_.num_threads;
#ifdef _OPENMP
    if (num_threads <= 0) num_threads = omp_get_max_threads();
    num_threads = std::min(num_threads, cfg_.tabu_attempts);
#else
    num_threads = 1;
#endif

    std::vector<std::vector<int>> results(cfg_.tabu_attempts);
    std::vector<bool> success(cfg_.tabu_attempts, false);

    // Deadline: leave some time for cleanup
    double remaining = cfg_.time_limit_seconds;
    {
        auto now = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(now - start_time_).count();
        remaining -= elapsed;
    }
    if (remaining < 0.1) return false;

    auto deadline = std::chrono::steady_clock::now() +
        std::chrono::duration_cast<std::chrono::steady_clock::duration>(
            std::chrono::duration<double>(remaining * 0.9));

    if (cfg_.verbose) {
        std::cout << "[XL] Launching " << cfg_.tabu_attempts << " tabu search workers"
                  << " on " << num_threads << " threads\n";
    }

#ifdef _OPENMP
    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
#endif
    for (int attempt = 0; attempt < cfg_.tabu_attempts; ++attempt) {
        results[attempt].resize(mn_);
        success[attempt] = tabu_worker(pair_sums, a, results[attempt],
                                       thread_seed + attempt * 7919,
                                       deadline);
    }

    // Find first successful attempt
    for (int attempt = 0; attempt < cfg_.tabu_attempts; ++attempt) {
        if (success[attempt]) {
            assignment_out = results[attempt];
            if (cfg_.verbose) {
                std::cout << "[XL] Worker " << attempt << " succeeded\n";
            }
            return true;
        }
    }
    return false;
}

// ============================================================
// Build result
// ============================================================
Result SolverXL::build_result(const std::vector<int>& M,
                               const std::vector<int>& N,
                               const std::vector<int>& assignment,
                               int a) const {
    Result r;
    r.found = true;
    r.a = a;
    r.set_m = M;
    r.set_n = N;
    int idx = 0;
    for (int u : M) {
        for (int v : N) {
            int label = assignment[idx];
            int weight = u + v + label;
            r.edge_labeling.push_back({{u, v}, label});
            r.edge_weights.push_back(weight);
            idx++;
        }
    }
    std::sort(r.edge_weights.begin(), r.edge_weights.end());
    return r;
}

// ============================================================
// Fall back to standard solver for small/medium graphs
// ============================================================
Result SolverXL::solve_with_standard_solver() {
    Config c;
    c.use_master_formula = cfg_.use_master_formula;
    c.use_ac3 = cfg_.use_ac3;
    c.use_greedy = cfg_.use_greedy;
    c.use_iterative_deepening = cfg_.use_iterative_deepening;
    c.use_mcv_ordering = cfg_.use_mcv_ordering;
    c.use_memoization = cfg_.use_memoization;
    c.use_restart = cfg_.use_restart;
    c.use_lcv_value_ordering = cfg_.use_lcv_value_ordering;
    c.use_target_propagation = cfg_.use_target_propagation;
    c.use_partition_ordering = cfg_.use_partition_ordering;
    c.use_local_search = cfg_.use_local_search;
    c.time_limit_seconds = cfg_.time_limit_seconds;
    c.per_partition_time = cfg_.per_partition_time;
    c.verbose = cfg_.verbose;

    Solver s(m_, n_, c);
    return s.solve();
}

// ============================================================
// Main solve entry point
// ============================================================
Result SolverXL::solve() {
    start_time_ = std::chrono::steady_clock::now();

    // Decide strategy: standard solver for small graphs, XL for large
    bool use_xl_mode = cfg_.skip_enumeration_xl && (mn_ > cfg_.enumeration_threshold);

    if (!use_xl_mode) {
        if (cfg_.verbose) {
            std::cout << "[XL] m*n=" << mn_ << " <= threshold " << cfg_.enumeration_threshold
                      << ", using standard 11-layer solver\n";
        }
        return solve_with_standard_solver();
    }

    if (cfg_.verbose) {
        std::cout << "[XL] m*n=" << mn_ << " > threshold, using parallel tabu search\n";
    }

    // Generate seed partitions
    std::vector<std::pair<std::vector<int>, std::vector<int>>> seeds;
    generate_seed_partitions(seeds);

    if (seeds.empty()) {
        if (cfg_.verbose) {
            std::cout << "[XL] No valid seed partitions for K_{" << m_ << "," << n_
                      << "}; falling back to standard 11-layer solver\n";
        }
        return solve_with_standard_solver();
    }

    if (cfg_.verbose) {
        std::cout << "[XL] Generated " << seeds.size() << " seed partitions\n";
    }

    Result result;
    int seed_idx = 0;
    for (auto& [M, N] : seeds) {
        if (time_exceeded()) {
            result.timed_out = true;
            result.status = "Time limit exceeded during partition search";
            break;
        }

        int a = compute_a(M, N);
        if (a < 0) continue;

        if (cfg_.verbose) {
            std::cout << "[XL] Trying seed partition #" << seed_idx
                      << " with a=" << a << "\n";
        }

        std::vector<int> pair_sums = compute_pair_sums(M, N);
        std::vector<int> assignment;

        bool ok = tabu_search(pair_sums, a, assignment, 0x9E3779B9 ^ seed_idx);
        if (ok) {
            result = build_result(M, N, assignment, a);
            result.status = "Solved by parallel tabu search (Layer 12)";
            result.partitions_tested = seed_idx + 1;
            break;
        }
        seed_idx++;
    }

    if (!result.found && !result.timed_out) {
        if (cfg_.verbose) {
            std::cout << "[XL] All " << seed_idx
                      << " seeds failed; falling back to standard 11-layer solver\n";
        }
        // Adjust remaining time budget
        auto now = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(now - start_time_).count();
        double remaining = cfg_.time_limit_seconds - elapsed;
        if (remaining > 1.0) {
            ConfigXL c2 = cfg_;
            c2.time_limit_seconds = remaining;
            SolverXL fallback(m_, n_, c2);
            // Direct call to standard solver
            Result fb = fallback.solve_with_standard_solver();
            if (fb.found) {
                fb.elapsed_ms = std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - start_time_).count();
                return fb;
            }
        }
        result.status = "No solution found (XL seeds + 11-layer fallback)";
    }

    auto end = std::chrono::steady_clock::now();
    result.elapsed_ms = std::chrono::duration<double, std::milli>(end - start_time_).count();
    return result;
}

} // namespace seatl
