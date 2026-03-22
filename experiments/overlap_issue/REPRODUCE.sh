#!/bin/bash
# ----------------------------------------------------------------
# Reprouduction script for HyTGraph data transfer overlap analysis
# ----------------------------------------------------------------

# 1. Generate large graph (1M nodes, 100M edges)
python3 gen_large_graph.py 1000000 100000000

# 2. Build the project
mkdir -p build && cd build
cmake ..
make -j 8

# 3. Profile with nsys (explicit mode, hybrid=1)
# Use 10 segments and 2 streams to see potential overlap
nsys profile -t cuda,nvtx -o bfs_profile_exp ./hybrid_bfs -graphfile ../large_graph.el -format market_big -SEGMENT 10 -n_stream 2 -hybrid 1 -source_node 0 --force-overwrite true

# 4. Extract GPU trace for analysis
nsys stats --report cuda_gpu_trace bfs_profile_exp.nsys-rep > gpu_trace.txt

echo "Experiment finished. See gpu_trace.txt for analysis."
