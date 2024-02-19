# reads2graph
Generate edit-distance-based read graph from short reads via looped minimiser and locality-sensitive hashing


## Dependencies
Sharg 1.1.1
Seqan 3.3.0
Boost 1.82.0
OpenMP 9.0.1
xxhash 0.8.2
omh 0.0.2 (source code of OMH compute used)

## Compile

```
cd build
cmake ..
make

../bin/reads2graph -i ../data/umi_SRR1543964.fastq -o ../result -k 15 -w 20 -x 6 -n 1 -g graph.dot
```
2024-01-11 17:53:22 - 2024-01-11 20:37:32 
../bin/reads2graph -i ../data/umi_SRR1543964.fastq -o ../result -k 15 -w 238 -x 6 -n 1 -g graph.dot
26 cores
Edit distance by minimiser-5: 3190363 pairs
Edit distance by minimiser-4: 2313668 pairs
Edit distance by minimiser-3: 1539273 pairs
Edit distance by minimiser-2: 1086861 pairs
Edit distance by minimiser-1: 102001 pairs

doubly using omp for pairwise comparisons
./bin/reads2graph -i ./data/umi_SRR1543964.fastq -o ./result -k 15 -w 238 -x 6 -n 1 -g graph2.dot