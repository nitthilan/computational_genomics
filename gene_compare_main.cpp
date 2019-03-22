//#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "fasta_read.h"
#include "gene_compare.h"






int main(int argc, char **argv){

	// 
	if(!(argc == 3 || argc == 4)){
		printf("<exe> <input file containing both s1 and s2> <0: global, 1: local> <optional: path to parameters config file>\n");
		printf("Num arguments %d\n", argc);
		exit(-1);
	}
	// Check input arguments
	char *p_filename = argv[1];
	int is_global = atoi(argv[2]);
	char *p_config_filename, config_filename[] = "parameters.config";
	if(argc == 4){
		p_config_filename = argv[3];
	}
	else{
		p_config_filename = config_filename;
	}

	// Check the input arguments
	if( access(p_filename, F_OK ) == -1 ) {
		printf("Input file: %s does not exist.\n", p_filename);
		exit(-1);
	} 
	if( access(p_config_filename, F_OK ) == -1 ) {
    	printf("Config file: %s does not exist.\n", p_config_filename);
		exit(-1);
	}
	if(!(is_global == 0 || is_global == 1)){
		printf("Wrong is_global flag %d\n", is_global);
		exit(-1);
	}

	// Read the input sequences
	// char *p_char_in_seq[2];
	input_data_t s_input_data;
	read_input_sequence(p_filename, &s_input_data);
	// Read config file
	param_config_t param_config;
	read_config(p_config_filename, &param_config);
	convert_to_small(&s_input_data);
	// printf("%s\n%s\n",s_input_data.p_char_in_seq[0], s_input_data.p_char_in_seq[1]);
	printf("Scores: match %d, mismatch %d, g %d, h %d\n", param_config.match, param_config.mismatch,
		param_config.g, param_config.h);
	printf("Sequence 1 = 's1', length = %d characters\n", s_input_data.num_char_in_seq[1]);
	printf("Sequence 2 = 's2', length = %d characters\n", s_input_data.num_char_in_seq[0]);

	dp_state_t dp_state;
	output_stat_t output_stat;
	find_alignment(&s_input_data, &dp_state, &param_config, is_global, &output_stat);
	print_output_stat(&output_stat);

	// free the allocated memory
	for(int i=0;i<2;i++){
		if(s_input_data.p_char_in_seq[i]){
			free(s_input_data.p_char_in_seq[i]);
		}
	}
	if(dp_state.p_dp_table){
		free(dp_state.p_dp_table);
	}
	if(output_stat.p_output_trace){
		free(output_stat.p_output_trace);
	}
}
