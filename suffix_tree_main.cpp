#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "fasta_read.h"
#include "suffix_tree.h"

// https://eecs.wsu.edu/~ananth/CptS571/Lectures/index.htm
// https://www.eecs.wsu.edu/~ananth/CptS571/Programs/Program2/index.htm

// ./suffix_tree_main ../data/suffix_tree/s1.fas ../data/suffix_tree/English_alphabet.txt
// ./suffix_tree_main ../data/suffix_tree/s2.fas ../data/suffix_tree/English_alphabet.tx


int main(int argc, char **argv){

	if(!(argc == 3)){
		printf("<exe> <input file containing s> <input alphabet file>\n");
		printf("Num arguments %d\n", argc);
		exit(-1);
	}
	// Check input arguments
	char *p_fast_filename = argv[1];
	char *p_alphabet_filename = argv[2];

	// Check the input arguments
	if( access(p_fast_filename, F_OK ) == -1 ) {
		printf("Fasta file: %s does not exist.\n", p_fast_filename);
		exit(-1);
	} 
	if( access(p_alphabet_filename, F_OK ) == -1 ) {
    	printf("Alphabet file: %s does not exist.\n", p_alphabet_filename);
		exit(-1);
	}

	// Read the input sequences
	// char *p_char_in_seq[2];
	input_data_t s_input_data;
	read_input_sequence(p_fast_filename, &s_input_data);
	// convert_to_small(&s_input_data);
	// Read alphabet file
	alphabets_t s_alphabets;
	// read_alphabets(p_alphabet_filename, &s_alphabets);

	// build a suffix tree
	st_node_t *p_st_root;
	st_construct(&s_input_data, &s_alphabets, &p_st_root);
	// st_construct_unoptimised(&s_input_data, &s_alphabets, &p_st_root);

	// free the allocated memory
	for(int i=0;i<2;i++){
		if(s_input_data.p_char_in_seq[i]){
			free(s_input_data.p_char_in_seq[i]);
		}
	}
}