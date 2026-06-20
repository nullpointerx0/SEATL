/*
 * seatl_cli_xl.cpp — CLI for the 12-layer XL SEATL solver.
 * Compatible with the standard backend parser.
 *
 * Usage:
 *   ./seatl_cli_xl M N [-t SECONDS] [--no-matrix]
 */

#include "seatl_solver_xl.hpp"
#include <iostream>
#include <iomanip>
#include <string>
#include <set>
#include <map>
#include <algorithm>
#include <vector>

using namespace seatl;

void print_matrix(const std::vector<int>& M, const std::vector<int>& N,
                  const std::vector<std::pair<std::pair<int,int>, int>>& edges) {
    std::map<std::pair<int,int>, int> weight_map;
    int total_weight = 0;
    for (auto& [uv, lbl] : edges) {
        int w = uv.first + uv.second + lbl;
        weight_map[uv] = w;
        total_weight += w;
    }
    std::cout << "\n  Edge Weight Matrix:\n";
    std::cout << "       ";
    for (int v : N) std::cout << std::setw(4) << v << " ";
    std::cout << "\n      ";
    for (size_t i = 0; i < N.size(); ++i) std::cout << "-----";
    std::cout << "\n";
    for (int u : M) {
        std::cout << "  " << std::setw(3) << u << " |";
        for (int v : N) {
            auto it = weight_map.find({u, v});
            if (it != weight_map.end()) std::cout << std::setw(4) << it->second << " ";
            else std::cout << "   - ";
        }
        std::cout << "\n";
    }
    std::cout << "\n  Total edge weight sum: " << total_weight << "\n";
}

int main(int argc, char** argv) {
    int m, n;
    double budget = 60.0;
    bool show_matrix = true;
    int threads = 0;

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " M N [-t SECONDS] [--no-matrix] [-j THREADS]\n";
        return 2;
    }
    try {
        m = std::stoi(argv[1]);
        n = std::stoi(argv[2]);
    } catch (...) {
        std::cerr << "m and n must be integers.\n";
        return 2;
    }

    for (int i = 3; i < argc; ++i) {
        std::string a = argv[i];
        if ((a == "-t" || a == "--time") && i + 1 < argc) budget = std::stod(argv[++i]);
        else if (a == "--no-matrix") show_matrix = false;
        else if (a == "--matrix") show_matrix = true;
        else if (a == "-j" && i + 1 < argc) threads = std::stoi(argv[++i]);
    }

    std::cout << "\n=== SEATL Solver (12 layers / XL) for K_{" << m << "," << n << "} ===\n";

    ConfigXL cfg;
    cfg.time_limit_seconds = budget;
    cfg.per_partition_time = budget;
    cfg.num_threads = threads;
    cfg.verbose = false;

    SolverXL s(m, n, cfg);
    Result r = s.solve();

    std::cout << "Time: " << std::fixed << std::setprecision(2) << r.elapsed_ms << " ms\n";
    std::cout << "Status: " << r.status << "\n";
    std::cout << "Diagnostics:\n";
    std::cout << "  Partitions tested: " << r.partitions_tested
              << " (L1 pruned: " << r.partitions_pruned_layer1
              << ", L2 pruned: " << r.partitions_pruned_layer2 << ")\n";
    std::cout << "  Layer 3 greedy successes: " << r.greedy_successes << "\n";
    std::cout << "  DFS invocations: " << r.dfs_invocations << "\n";
    std::cout << "  Memo hits: " << r.memo_hits << "\n";

    if (!r.found) {
        std::cout << "\nNo solution found.\n";
        return 1;
    }

    std::cout << "\n--- Solution (a = " << r.a
              << ", AP = " << r.a << ".." << (r.a + m*n - 1) << ") ---\n";
    std::cout << "  M = {";
    for (size_t i = 0; i < r.set_m.size(); ++i)
        std::cout << r.set_m[i] << (i + 1 < r.set_m.size() ? "," : "");
    std::cout << "}    N = {";
    for (size_t i = 0; i < r.set_n.size(); ++i)
        std::cout << r.set_n[i] << (i + 1 < r.set_n.size() ? "," : "");
    std::cout << "}\n";

    if (show_matrix) print_matrix(r.set_m, r.set_n, r.edge_labeling);

    std::cout << "  All edges (" << r.edge_labeling.size() << " total):\n";
    for (auto& [uv, lbl] : r.edge_labeling) {
        int w = uv.first + uv.second + lbl;
        std::cout << "    (" << std::setw(3) << uv.first
                  << "," << std::setw(3) << uv.second << ")  label="
                  << std::setw(4) << lbl
                  << "  weight=" << std::setw(4) << w << "\n";
    }

    return 0;
}
