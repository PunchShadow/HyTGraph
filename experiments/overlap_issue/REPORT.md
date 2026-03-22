# HyTGraph Data Transfer Overhead and Overlap Analysis

## 1. Analysis Goal
Analyze the data transfer overhead in HyTGraph, particularly the ratio of active nodes to overhead, and determine if data transfer can be overlapped with kernel computation.

## 2. Methodology
- **Datasets:** 
  1. Synthetic random graph: 1,000,000 nodes, 100,000,000 edges.
  2. **Orkut Simulation:** 3,000,000 nodes, 117,000,000 edges.
  3. **uk-2007 (Real World):** 105,153,953 nodes, 3,301,876,564 edges (~13GB).
- **Algorithm:** BFS (from `samples/hybrid_bfs`).
- **Configuration:** 
  - `SEGMENT=10` or `20` (splitting the graph to manage memory).
  - `n_stream=2` (use 2 CUDA streams for concurrent execution).
  - `hybrid=1` (use explicit data transfer for all segments).
- **Tooling:** `nsys` (NVIDIA Nsight Systems) for capturing execution timeline and GPU traces.

## 3. Experimental Findings

### 3.1 Data Transfer Overhead
In **Explicit Mode** (`hybrid=1`), the system transfers the entire segment (all its edges) from host to device if *at least one* node in that segment is active.
- **Orkut Scale (117M edges):** Each segment (1/10) contains ~11.7M edges (46.8MB). Single segment transfer takes **~3.9ms**.
- **uk-2007 Scale (3.3B edges):** With `SEGMENT=20`, each segment contains ~165,000,000 edges.
  - Each edge is 4 bytes, resulting in a **~660MB** transfer per segment per active round.
  - At PCIe Gen3 x16 (~12GB/s), a single segment transfer takes **~55ms**.
  - **Case Study (uk-2007):** In a 53-round BFS, we observed 1,040 segment activations. 
    - Theoretical transfer time: 1040 * 55ms = **57.2 seconds**.
    - Observed Total Time: **63.4 seconds**.
    - Data transfer overhead accounts for **~90%** of total execution time for massive graphs.

### 3.2 Active Nodes/Frontiers to Overhead Ratio
- **Explicit Mode:** The transfer size is constant (full segment size) regardless of frontier density.
  - **Inefficiency:** In `uk-2007`, Round 1 with `source_node=1` (only 16 neighbors) still triggered a full **660MB** transfer.
  - **Efficiency:** The ratio is best during the "peak" of BFS when frontiers are dense across segments.
- **Compaction Mode:** (Based on code analysis) Only active edges are transferred, making overhead proportional to the frontier size.

### 3.3 Overlap of Data Transfer and Computation
Using `nsys` profile, we confirmed:
1. **Concurrency:** Multiple CUDA streams are utilized effectively.
2. **Overlap Proof:** 
   - While one stream is busy with a 3.9ms - 55ms `memcpy` (H2D), another stream executes the `SyncPushDDB` kernel.
   - For BFS, the computation is so fast (tens of microseconds) that it is effectively **hidden** by the multi-millisecond data transfers.
3. **Efficiency:** `n_stream > 1` is critical. It allows the PCIe bus to stay saturated by pre-fetching the next segment's data while the GPU is processing the current one.

## 4. Conclusion
HyTGraph successfully implements **data transfer and computation overlap**. However, performance on large-scale graphs is strictly **bandwidth-bound**. The "active node to overhead ratio" in Explicit mode is inversely proportional to frontier density, leading to high overhead in the early/late stages of graph traversal.

## 5. Reproduction
Run the following script from the project root:
```bash
bash experiments/overlap_issue/REPRODUCE.sh
```
