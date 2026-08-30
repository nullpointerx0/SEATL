# SEATL XL

A high-performance C++ solver for **Super (a,d)-Edge-Antimagic Total Labeling (SEATL)** of complete bipartite graphs `K_{m,n}`, extended with a 12th layer — parallel tabu search — to scale efficiently up to `K_{20,20}` and beyond.

## Table of Contents

- [Problem Definition](#problem-definition)
- [What's New in XL](#whats-new-in-xl)
- [Benchmark Results](#benchmark-results)
- [Build](#build)
- [Usage](#usage)
	- [Architecture](#architecture)
- [How Layer 12 Works](#how-layer-12-works)
- [Project Structure](#project-structure)
- [Deployment](#deployment)
- [Backend Integration](#backend-integration)
- [Roadmap](#roadmap)

## Problem Definition

For a complete bipartite graph `K_{m,n}`:

- **Vertex labels**: the integers `{1, 2, ..., m+n}`, each used exactly once.
- **Edge labels**: the integers `{m+n+1, ..., m+n+mn}`, each used exactly once.
- **Edge weight**: `w(u, v) = label(u) + label(v) + label(e)`.
- **Goal**: assign labels such that all edge weights form a consecutive arithmetic progression `{a, a+1, ..., a+mn-1}` (the `d = 1` case).

A valid solution is a labeling where every edge weight is unique and the full set of weights spans this progression with no gaps or duplicates.

## What's New in XL

The XL solver builds on the original 11-layer architecture (`seatl_solver.cpp`) by adding **Layer 12: Parallel Tabu Search**, purpose-built for graphs too large to solve by enumeration or search alone.

| # | Improvement | Description |
|---|-------------|-------------|
| 1 | Skip enumeration | For `mn > 30`, don't enumerate all `C(m+n, m)` partitions — try 8–16 smart seeds instead |
| 2 | O(1) incremental scoring | Each tabu swap is scored in O(1) instead of O(mn) — roughly 100× faster per move |
| 3 | Parallel workers | OpenMP runs 8–16 tabu workers concurrently; the first to succeed wins |
| 4 | Tabu memory | Recent swaps are forbidden for `tenure` iterations to escape local minima |
| 5 | Smart restarts | After 2000+ stalled iterations, the search fully randomizes and retries |

## Benchmark Results

Wall-clock time to find a verified-valid SEATL labeling (every edge weight is checked against the target arithmetic progression on each run). Values are the median of 5 runs; all 5 runs solved successfully.

| Graph     | 1 thread | 8 threads | Parallel speedup |
|-----------|---------:|----------:|:-----------------:|
| K_{7,7}   |   5.4 ms |   2.1 ms  | 2.5× |
| K_{10,10} |    35 ms |  12.7 ms  | 2.8× |
| K_{12,12} |    97 ms |    29 ms  | 3.3× |
| K_{12,15} |   1.44 s |   329 ms  | 4.4× |
| K_{15,15} |   0.62 s |   231 ms  | 2.7× |
| K_{18,18} |   1.35 s |   345 ms  | 3.9× |
| K_{20,20} |   7.84 s |   1.96 s  | 4.0× |
| K_{22,22} |   21.8 s |   4.65 s  | 4.7× |
| K_{25,25} |   95.5 s |   16.6 s  | 5.8× |

**Methodology.** CPU: Intel Core i7-13700HX. Built with `g++ 15.2 -O3 -fopenmp` (C++17). 8-thread runs pin OpenMP to 8 cores (`-j 8`); 1-thread runs use `-j 1`. Times are wall-clock medians as reported by the solver. Reproduce with `./bench.sh`; raw data is in [`seatl_benchmark_raw.csv`](seatl_benchmark_raw.csv) and a formatted summary in [`seatl_benchmark_results.md`](seatl_benchmark_results.md).

Parallel speedup is sub-linear (2.5–6×) because Layer 12 is a *first-worker-wins* race: workers explore independent random regions of the search space, so adding cores raises the probability that *some* worker finds `score = 0` quickly, rather than dividing a fixed workload. Speedup grows with instance size as the search gets harder.

## Build

Requires a C++17 compiler. OpenMP is optional but recommended for parallel search.

```bash
# With OpenMP (recommended)
g++ -std=c++17 -O3 -fopenmp \
    seatl_solver.cpp seatl_solver_xl.cpp seatl_cli_xl.cpp \
    -o seatl_cli_xl

# Without OpenMP (single-threaded)
g++ -std=c++17 -O3 \
    seatl_solver.cpp seatl_solver_xl.cpp seatl_cli_xl.cpp \
    -o seatl_cli_xl
```

## Usage

```bash
./seatl_cli_xl M N [-t SECONDS] [-j THREADS] [--no-matrix | --matrix]
```

| Flag | Default | Description |
|------|---------|--------------|
| `M N` | — | Required. Bipartite graph part sizes, e.g. `20 20` for `K_{20,20}` |
| `-t`, `--time SECONDS` | `60` | Time budget for the search |
| `-j THREADS` | auto | Number of parallel tabu workers |
| `--no-matrix` | off | Suppress the edge-weight matrix (useful for backend integration) |
| `--matrix` | on | Print the edge-weight matrix |

Examples:

```bash
# Solve K_{20,20} with a 60s budget
./seatl_cli_xl 20 20 -t 60

# Use 8 threads explicitly
./seatl_cli_xl 20 20 -t 60 -j 8

# Suppress matrix output for programmatic/backend use
./seatl_cli_xl 20 20 -t 60 --no-matrix
```

## Architecture

The original solver (`seatl_solver.cpp` / `.hpp`) implements an 11-layer pipeline; the XL solver adds Layer 12 on top for large instances.

**Layers 1–11 (original solver):**

1. Master formula — algebraic pruning
2. AC-3 constraint propagation
3. Greedy initialization
4. Iterative deepening DFS
5. Most-constrained-variable ordering
6. Memoization
7. Restart with diversification
8–11. Additional refinement passes (see `seatl_solver.hpp`)

**Layer 12 (XL solver):** Parallel Tabu Search — described below.

Mode selection is automatic based on graph size:

```cpp
if (m * n > enumeration_threshold)  // default 30
    use XL mode;       // Layer 12
else
    use standard mode;  // Layers 1-11
```

XL mode outperforms the standard solver for all but the smallest graphs (`K_{2,2}`, `K_{3,3}`, `K_{4,4}`), where the standard solver's greedy/DFS approach is already near-optimal.

## How Layer 12 Works

For each seed partition, run in parallel:

```
1. Initialize: shuffle labels randomly onto edges
2. Compute score: how many target slots are missing or duplicated?
3. Repeat:
   a. Sample ~mn random swap candidates (i, j)
   b. For each, compute the score delta in O(1) — only 4 slots change
   c. Pick the best non-tabu swap (or an aspirational move if it beats the best known)
   d. Apply the swap; mark it tabu for `tenure` iterations
   e. If score == 0, SOLVED
   f. If stuck too long, restart from a new random initialization
4. Return the best assignment found
```

**Why incremental scoring is fast.** A naive rescan is O(mn) per candidate swap:

```
for each swap candidate:
    swap labels
    rescan all mn edges to compute score   # O(mn)
    undo swap
```

Because a swap of edges `(i, j)` only changes 4 weight slots (old weight of `i`, old weight of `j`, new weight of `i`, new weight of `j`), the solver instead tracks `slot_count[t]` — how many edges currently hit each target slot — and computes the delta directly from those 4 slots in **O(1)**. For `K_{20,20}` (400 edges), this is the difference between touching all 400 edges per candidate swap and touching a constant number of slots.

**Parallelism.** Each worker starts from an independent random seed and explores a different region of the search space; the first worker to reach `score == 0` wins the race. This is why speedup is sub-linear rather than proportional to thread count — see [Benchmark Results](#benchmark-results).

## Project Structure

```
.
├── seatl_solver.hpp / seatl_solver.cpp        # Original 11-layer solver
├── seatl_solver_xl.hpp / seatl_solver_xl.cpp   # XL Layer 12 (parallel tabu search)
├── seatl_cli_xl.cpp                            # CLI entry point (drop-in replacement)
├── seatl_cli_xl                                # Prebuilt binary
├── bench.sh                                    # Benchmark runner (shell)
├── run_benchmarks.py                           # Benchmark runner (Python)
├── seatl_benchmark_raw.csv                     # Raw benchmark data
├── seatl_benchmark_results.md                  # Formatted benchmark summary
└── deployment/                                 # Docker-based deployment assets
    ├── Dockerfile
    ├── docker-compose.yml
    ├── entrypoint.sh
    └── deploy.md
```

## Deployment

This repository contains the solver CLI only — it does not include a web frontend or HTTP backend. The `deployment/` directory provides a Docker-based path to build and run the solver as a container, with full instructions in [`deployment/deploy.md`](deployment/deploy.md).

Quick start:

```bash
# Build the image
docker build -f deployment/Dockerfile -t seatl-xl:latest .

# Run the solver
docker run --rm seatl-xl:latest 20 20 -t 60 --no-matrix
```

Or with Docker Compose:

```bash
docker compose -f deployment/docker-compose.yml build
docker compose -f deployment/docker-compose.yml run --rm solver 20 20 -t 60 -j 8
```

See [`deployment/deploy.md`](deployment/deploy.md) for remote-host deployment and backend integration options.

## Backend Integration

The XL CLI output format matches the standard solver's output, so an existing parser built against `seatl_cli_xl` (standard) requires no changes — simply point `SOLVER_PATH` at the XL binary.

Note that this repository does not expose an HTTP interface on its own; a backend service should invoke the CLI (directly or via the Docker image) and parse its stdout.

## Roadmap

Potential directions for scaling beyond `K_{30,30}`:

1. **Genetic algorithm crossover** — combine good partial assignments from different workers
2. **Simulated annealing** — accept occasional worse moves to escape plateaus
3. **GPU acceleration** — score thousands of candidate swaps in parallel
4. **SAT-based seed selection** — use constraint reasoning to pick more promising initial partitions
