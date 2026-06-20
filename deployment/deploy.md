# Deployment Guide

This repository currently contains a **C++ solver CLI**, not a complete web app. The README mentions a Node.js backend parser, but that backend is not present in this repository.

The files in this `deployment/` folder let you:

1. Build the solver into a Docker image.
2. Run it locally or on a remote Linux host.
3. Reuse the container from another backend by invoking the CLI inside the container.

If you intended to deploy a web app, you will also need the missing backend/frontend repository.

## Included Files

- `Dockerfile`: multi-stage build for `seatl_cli_xl`
- `docker-compose.yml`: simple service definition for local or remote Docker Compose use
- `entrypoint.sh`: usage wrapper for the solver container

## Prerequisites

On the remote host, install:

- Docker Engine
- Docker Compose plugin (`docker compose`)

Verify:

```bash
docker --version
docker compose version
```

## Build the Image

From the repository root:

```bash
docker build -f deployment/Dockerfile -t seatl-xl:latest .
```

Or with Compose:

```bash
docker compose -f deployment/docker-compose.yml build
```

## Run the Solver

Direct Docker:

```bash
docker run --rm seatl-xl:latest 20 20 -t 60 --no-matrix
```

With Compose:

```bash
docker compose -f deployment/docker-compose.yml run --rm solver 20 20 -t 60 --no-matrix
```

Use explicit threads:

```bash
docker compose -f deployment/docker-compose.yml run --rm solver 20 20 -t 60 -j 8
```

## Deploy on a Remote Host

### Option 1: Copy the repository to the server

On your local machine:

```bash
scp -r /path/to/EL_Crazy user@your-server:/opt/EL_Crazy
```

On the server:

```bash
cd /opt/EL_Crazy
docker compose -f deployment/docker-compose.yml build
docker compose -f deployment/docker-compose.yml run --rm solver 20 20 -t 60 --no-matrix
```

### Option 2: Pull from git on the server

On the server:

```bash
git clone <your-repository-url> /opt/EL_Crazy
cd /opt/EL_Crazy
docker compose -f deployment/docker-compose.yml build
docker compose -f deployment/docker-compose.yml run --rm solver 20 20 -t 60 --no-matrix
```

## Integrating With a Backend

If your actual web app lives in another repository, point that backend to this containerized solver in one of these ways:

1. Build and run this image separately, then have the backend call `docker run --rm seatl-xl:latest ...`.
2. Copy the same `Dockerfile` into the backend repo and build the solver into the backend image.
3. Use a multi-service Compose setup where the backend service shells out to this image.

Because this repo does not expose HTTP by itself, it is **not** a standalone web service container.

## Recommended Production Follow-Up

If you want a real remote-host web deployment, the next step is to add one of these:

1. A small HTTP wrapper around `seatl_cli_xl` in Node.js, Go, or Python.
2. A combined `docker-compose.yml` with the backend, reverse proxy, and this solver image.
3. Environment-specific logging, timeouts, and request validation.
