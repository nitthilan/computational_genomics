cmake_minimum_required (VERSION 2.6)
project (computational_genomics)

# add the executable
add_executable (suffix_tree_main suffix_tree_main.cpp suffix_tree.cpp suffix_tree.h fasta_read.cpp fasta_read.h)
add_executable (gene_compare_main gene_compare_main.cpp gene_compare.cpp gene_compare.h fasta_read.cpp fasta_read.h)

add_executable (read_mapping_main_1 read_mapping_main.cpp read_mapping.cpp read_mapping.h suffix_tree.cpp suffix_tree.h gene_compare.cpp gene_compare.h fasta_read.cpp fasta_read.h)
