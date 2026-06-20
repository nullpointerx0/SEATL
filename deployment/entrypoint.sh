#!/bin/sh
set -eu

if [ "$#" -eq 0 ] || [ "$1" = "--help" ]; then
  cat <<'EOF'
Usage:
  docker run --rm seatl-xl M N [-t SECONDS] [--no-matrix] [-j THREADS]

Examples:
  docker run --rm seatl-xl 20 20 -t 60 --no-matrix
  docker run --rm seatl-xl 20 20 -t 60 -j 8
EOF
  exit 0
fi

if [ -n "${SEATL_DEFAULT_ARGS}" ]; then
  # shellcheck disable=SC2086
  exec /usr/local/bin/seatl_cli_xl ${SEATL_DEFAULT_ARGS} "$@"
fi

exec /usr/local/bin/seatl_cli_xl "$@"
