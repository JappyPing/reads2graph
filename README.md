# ReadGraph
Generate edit-distance-based read graph from short reads


## Complie

```
cd build
cmake ..
make

../bin/ReadGraph -i ../data/umi_SRR1543964.fastq -o ../result -k 15 -w 20 -x 6 -n 1 -g graph.dot
```
2024-01-11 17:53:22 - 2024-01-11 20:37:32 
../bin/ReadGraph -i ../data/umi_SRR1543964.fastq -o ../result -k 15 -w 238 -x 6 -n 1 -g graph.dot
26 cores
Edit distance by minimiser-5: 3190363 pairs
Edit distance by minimiser-4: 2313668 pairs
Edit distance by minimiser-3: 1539273 pairs
Edit distance by minimiser-2: 1086861 pairs
Edit distance by minimiser-1: 102001 pairs

doubly using omp for pairwise comparisons
./bin/ReadGraph -i ./data/umi_SRR1543964.fastq -o ./result -k 15 -w 238 -x 6 -n 1 -g graph2.dot