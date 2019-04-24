//https://eecs.wsu.edu/~ananth/CptS571/Programs/Program3/index.htm



#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "fasta_read.h"
#include "suffix_tree.h"
#include "read_mapping.h"

// https://eecs.wsu.edu/~ananth/CptS571/Lectures/index.htm
// https://www.eecs.wsu.edu/~ananth/CptS571/Programs/Program2/index.htm

// ./read_mapping_main ../data/hw3/Peach_reference.fasta ../data/hw3/Peach_simulated_reads.fasta
// ./read_mapping_main ../data/hw3/reference.fasta ../data/hw3/read.fasta

int main(int argc, char **argv){

	if(!(argc == 3)){
		printf("<exe> <input file containing genome sequence G> <input file containing reads>\n");
		printf("Num arguments %d\n", argc);
		exit(-1);
	}
	// Check input arguments
	char *p_fast_filename = argv[1];
	char *p_read_filename = argv[2];

	// Check the input arguments
	if( access(p_fast_filename, F_OK ) == -1 ) {
		printf("Genome Sequence Fasta file: %s does not exist.\n", p_fast_filename);
		exit(-1);
	} 
	if( access(p_read_filename, F_OK ) == -1 ) {
    	printf("Read Fasta file: %s does not exist.\n", p_read_filename);
		exit(-1);
	}

	read_mapping(p_fast_filename, p_read_filename);
}