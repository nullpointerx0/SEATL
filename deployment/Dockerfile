FROM gcc:14-bookworm AS builder

WORKDIR /app

COPY seatl_solver.cpp seatl_solver.hpp ./
COPY seatl_solver_xl.cpp seatl_solver_xl.hpp ./
COPY seatl_cli_xl.cpp ./

RUN g++ -std=c++17 -O3 -fopenmp \
    seatl_solver.cpp seatl_solver_xl.cpp seatl_cli_xl.cpp \
    -o /out/seatl_cli_xl

FROM debian:bookworm-slim

RUN apt-get update \
    && apt-get install -y --no-install-recommends libgomp1 ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY --from=builder /out/seatl_cli_xl /usr/local/bin/seatl_cli_xl
COPY deployment/entrypoint.sh /entrypoint.sh

RUN chmod +x /entrypoint.sh

ENV SEATL_DEFAULT_ARGS=""

ENTRYPOINT ["/entrypoint.sh"]
CMD ["--help"]
