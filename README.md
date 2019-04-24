# computational_genomics
Project as part of assignment for computational genomics

```
mkdir build
cmake ../
make

```
- Produces three exe "gene_compare_main", "suffix_tree_main", "read_mapping_main_1"
- The exe "read_mapping_main_1" can be run by 

```
./read_mapping_main_1 ../data/hw3/Peach_reference.fasta ../data/hw3/Peach_simulated_reads.fasta 
```

Folder paths:
- build - Build and exe generated for testing
- data - All the data files used for testing

Files:
- fast_read.cpp, fasta_read.h - files for reading fasta files
- gene_compare.cpp/.h - global and local alignment function
- suffix_tree.cpp/.h - suffix tree generation and related search functions
- read_mapping.cpp/.h - read_mapping related files