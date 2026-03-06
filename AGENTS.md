# Repository Guidelines

## Project Structure & Module Organization
- `CMakeLists.txt`: top-level build configuration (CUDA + C++).
- `src/`: core implementation (parsers, utilities, etc.).
- `include/`: public headers and framework code.
- `samples/`: runnable apps (`hybrid_bfs`, `hybrid_cc`, `hybrid_pr`, `hybrid_sssp`) and small tools.
- `deps/`: third-party dependencies (e.g., gflags, cub, json).
- `datasets/SNAP/`: SNAP input graphs and helpers (optional).
- `build/`: out-of-source build directory.

## Build, Test, and Development Commands
- Configure + build:
  - `mkdir -p build && cd build`
  - `cmake .. -DCMAKE_BUILD_TYPE=Release`
  - `make -j`
- Run examples:
  - `./hybrid_bfs -graphfile ../example.el -format market_big -weight_num 1`
  - `./hybrid_bfs -graphfile ../datasets/SNAP/as-skitter.txt -format snap -weight_num 1`
- Help (gflags):
  - `./hybrid_bfs --help` (note: `-h` is not supported)

## Coding Style & Naming Conventions
- Language: CUDA C++ (C++11).
- Indentation: 4 spaces; braces on their own line are common in existing files.
- Naming:
  - Types/classes: `PascalCase` (e.g., `HybridEngine`).
  - Methods/variables: `camelCase` or `snake_case` as already used in the file.
  - Flags: lowercase with underscores (gflags), e.g., `--source_node`.
- No formatter or linter is configured; follow nearby file style.

## Testing Guidelines
- No dedicated unit-test framework is present.
- Use app-level checks when available:
  - `--check` compares against a CPU baseline (may be slow).
- Smoke test via running small graphs from `example.el` or `datasets/SNAP/`.

## Commit & Pull Request Guidelines
- Recent history uses simple, imperative messages like `Update README.md` or `Create example.el`.
- Recommendation: keep messages short and specific to the change.
- PRs should include:
  - Summary of changes and rationale.
  - Repro commands (build/run).
  - Any dataset/format assumptions (e.g., `-format snap`, `-undirected`).

## Data & Formats
- Supported formats: `market_big` and `snap`.
- For unweighted graphs, use `-weight_num 1` to generate weights.
- For undirected SNAP datasets (`*.ungraph.txt`), pass `--undirected`.
