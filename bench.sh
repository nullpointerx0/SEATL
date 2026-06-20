#!/usr/bin/env bash
# Benchmark + correctness harness for SEATL XL solver.
# For each K_{m,n}: run TRIALS times, verify weights == {a..a+mn-1}, report median ms.
set -u
BIN=./seatl_bench.exe
TRIALS=5
BUDGET=120

sizes=("7 7" "10 10" "12 12" "12 15" "15 15" "18 18" "20 20" "22 22" "25 25")

printf "%-10s | %-8s | %-10s | %-10s | %-8s\n" "Graph" "Verify" "median_ms" "min_ms" "solved/N"
printf -- "-----------+----------+------------+------------+---------\n"

for sz in "${sizes[@]}"; do
  m=$(echo $sz | cut -d' ' -f1); n=$(echo $sz | cut -d' ' -f2)
  mn=$((m*n))
  times=(); solved=0; verify="OK"
  for ((t=0; t<TRIALS; t++)); do
    out=$($BIN $m $n -t $BUDGET --no-matrix 2>/dev/null)
    ms=$(echo "$out" | grep -oP 'Time:\s*\K[0-9.]+')
    if echo "$out" | grep -q "All edges"; then
      solved=$((solved+1))
      a=$(echo "$out" | grep -oP 'a = \K[0-9]+' | head -1)
      # extract weights, verify they equal {a..a+mn-1}
      ok=$(echo "$out" | grep -oP 'weight=\s*\K[0-9]+' | sort -n | uniq | wc -l)
      lo=$(echo "$out" | grep -oP 'weight=\s*\K[0-9]+' | sort -n | head -1)
      hi=$(echo "$out" | grep -oP 'weight=\s*\K[0-9]+' | sort -n | tail -1)
      if [ "$ok" -ne "$mn" ] || [ "$lo" != "$a" ] || [ "$hi" != "$((a+mn-1))" ]; then
        verify="FAIL"
      fi
      times+=("$ms")
    fi
  done
  if [ ${#times[@]} -eq 0 ]; then
    printf "%-10s | %-8s | %-10s | %-10s | %d/%d\n" "K_{$m,$n}" "NONE" "-" "-" "$solved" "$TRIALS"
    continue
  fi
  sorted=($(printf '%s\n' "${times[@]}" | sort -n))
  cnt=${#sorted[@]}
  med=${sorted[$((cnt/2))]}
  mn_t=${sorted[0]}
  printf "%-10s | %-8s | %-10s | %-10s | %d/%d\n" "K_{$m,$n}" "$verify" "$med" "$mn_t" "$solved" "$TRIALS"
done
