
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>


#include "fasta_read.h"
#include "gene_compare.h"
#include "suffix_tree.h"
#include "read_mapping.h"

void input_data_alloc(input_data_t *p_input_data, int num_read_seq){
	//
	p_input_data->pp_char_in_seq = (char **)malloc(num_read_seq*sizeof(char*));
	p_input_data->p_num_char_in_seq = (int *)malloc(num_read_seq*sizeof(int));
	p_input_data->p_num_lines_in_seq = (int *)malloc(num_read_seq*sizeof(int));

	if(p_input_data->pp_char_in_seq == 0 || p_input_data->p_num_char_in_seq == 0 
	|| p_input_data->p_num_lines_in_seq == 0){
		printf("Memory allocation for input data failed. %d \n",
			num_read_seq);
		exit(-1);
	}
}
void input_data_dealloc(input_data_t *p_input_data, int num_read_seq){
	for(int i=0;i<num_read_seq;i++){
		if(p_input_data->pp_char_in_seq[i]){
			free(p_input_data->pp_char_in_seq[i]);
		}
	}
	if(p_input_data->pp_char_in_seq){
		free(p_input_data->pp_char_in_seq);
	}
	if(p_input_data->p_num_lines_in_seq){
		free(p_input_data->p_num_lines_in_seq);
	}
	if(p_input_data->p_num_char_in_seq){
		free(p_input_data->p_num_char_in_seq);
	}
}


// ./read_mapping_main ../data/hw3/Peach_reference.fasta ../data/hw3/Peach_simulated_reads.fasta
// ./read_mapping_main ../data/hw3/reference.fasta ../data/hw3/read.fasta

void read_mapping(char *p_reference_fasta, char *p_read_fasta){
	// Read the input sequences
	input_data_t s_input_data;
	input_data_alloc(&s_input_data, 1);
	input_data_t s_compare_data;
	input_data_alloc(&s_compare_data, 2);


	read_input_sequence(p_reference_fasta, &s_input_data, 0, 1);
	printf("Reference: number of characters %d, %d\n", s_input_data.p_num_char_in_seq[0],
		s_input_data.p_num_lines_in_seq[0]);
	// build a suffix tree
	// Read alphabet file
	alphabets_t s_alphabets;
	// read_alphabets(p_alphabet_filename, &s_alphabets);	
	st_node_t *p_st_root;
	st_construct(&s_input_data, &s_alphabets, &p_st_root);

	// Modify the suffix tree to store the leaf index
	int *p_leaf_map_id = (int *)malloc(s_input_data.p_num_char_in_seq[0]*sizeof(int));
	if(p_leaf_map_id == 0){
		printf("Leaf map index memory allocation failed. %d \n",s_input_data.p_num_char_in_seq[0]);
		exit(-1);
	}
	int cur_arr_idx = 0;

	clock_t begin_prep_st = clock();
	init_leaf_idx_range(p_st_root, p_leaf_map_id, &cur_arr_idx);
	clock_t end_prep_st = clock();
	double time_spent_st = (double)(end_prep_st - begin_prep_st) / CLOCKS_PER_SEC;
	printf("Time consumed prepare st %f\n", time_spent_st);

	// int max_depth;
	// dfs(p_st_root, &s_input_data, &max_depth);

	//
	// printf("ST leaves ");
	// for(int i=0; i<s_input_data.p_num_char_in_seq[0]; i++){
	// 	printf("%d, ", p_leaf_map_id[i]);
	// }
	// printf("\n");

	// Get the total number of sequences in reading file
	int num_read_seq = get_num_fasta(p_read_fasta);
	printf("Read: num_seq %d, %s\n", num_read_seq, p_read_fasta);

	// Read all the input strings to be checked
	input_data_t s_read_fasta;
	input_data_alloc(&s_read_fasta, num_read_seq);

	// input_data_t s_read_fasta;
	// s_read_fasta.pp_char_in_seq = (char **)malloc(num_read_seq*sizeof(char*));
	// s_read_fasta.p_num_char_in_seq = (int *)malloc(num_read_seq*sizeof(int));
	// s_read_fasta.p_num_lines_in_seq = (int *)malloc(num_read_seq*sizeof(int));
	// if(s_read_fasta.pp_char_in_seq == 0 || s_read_fasta.p_num_char_in_seq == 0 
	// || s_read_fasta.p_num_lines_in_seq == 0){
	// 	printf("Memory allocation for Read Fasta String failed. %d \n",
	// 		num_read_seq);
	// 	exit(-1);
	// }

	// Read all the input sequences
	read_input_sequence(p_read_fasta, &s_read_fasta, 0, num_read_seq);
	// for(int i=0; i<num_read_seq; i++){
	// 	printf("Size of the read fasta %d, %d\n", i, s_read_fasta.p_num_char_in_seq[i]);
	// }

	int num_local_search_options = 0;
	int num_misses = 0;
	clock_t begin = clock();
	for(int i=0; i<num_read_seq; i++){
		// get the corresponding read sequence
		// check file read of sequences require memory freeing
		// input_data_t s_read_fasta;
		// read_input_sequence(p_read_fasta, &s_read_fasta, i, 1);
		// printf("Size of the read sequence %d, %d\n", i, s_read_fasta.p_num_char_in_seq[i]);

		// if(i%1000 == 1){
		// 	printf("Num sequence is %d\n", i);
		// }

		// find the exact matching subsequence for the possible read sequence
		// check how to add matching for 25 continous characters
		// Note: Probably return a list of matching nodes
		int min_number_of_matches = 25;
		// int min_number_of_matches = 4;
		st_node_t *p_matching_node = find_node_matching_path(p_st_root, 
			s_input_data.pp_char_in_seq[0], s_read_fasta.pp_char_in_seq[i], 
			0, s_read_fasta.p_num_char_in_seq[i], min_number_of_matches);

		// if(p_matching_node->string_depth <= 25){
		// 	printf("MAtching unecessary node %d\n", p_matching_node->string_depth);
		// }


		// st_node_t *p_matching_node_sl = find_node_matching_path_sl(p_st_root, 
		// 	s_input_data.pp_char_in_seq[0], s_read_fasta.pp_char_in_seq[i], 
		// 	0, s_read_fasta.p_num_char_in_seq[i], min_number_of_matches);

		// if(p_matching_node_sl != p_matching_node){
		// 	printf("Error in suffix link %d, %d\n", p_matching_node_sl->node_idx, p_matching_node->node_idx);
		// }

		// if(p_matching_node){
		// 	printf("Matching suffixes ");
		// 	for (int j = p_matching_node->leaf_idx_range[0]; 
		// 			j <= p_matching_node->leaf_idx_range[1]; j++){
		// 		printf("%d, ", p_leaf_map_id[j]);
		// 	}
		// 	printf("\n");
		// }

		if(!p_matching_node){
			// printf("Error in matching %d\n", i);
			printf("READ_%d_from_Peach_reference.fasta No hit found\n", i+1);
			num_misses += 1;
			continue;
		}
		num_local_search_options += (p_matching_node->leaf_idx_range[1] - p_matching_node->leaf_idx_range[0]);


		int best_offset = 0;
		float max_length_coverage = 0;
		float max_percent_identity = 0;
		int best_align_length = 0;
		// for each possible matches calculate local matching score
		for(int j=p_matching_node->leaf_idx_range[0]; j<=p_matching_node->leaf_idx_range[1]; j++){

			int start_leaf_idx = p_leaf_map_id[j] - s_read_fasta.p_num_char_in_seq[i];
			int total_sequence_length = 2*s_read_fasta.p_num_char_in_seq[i];
			if(start_leaf_idx<0){
				start_leaf_idx = 0;
			}
			if(start_leaf_idx + total_sequence_length > s_input_data.p_num_char_in_seq[0]){
				total_sequence_length = s_input_data.p_num_char_in_seq[0] - start_leaf_idx;
			}

			if(total_sequence_length <= 50){
				printf("Sequence length small %d, %d, %d\n", i, p_leaf_map_id[j], total_sequence_length);
			}

			// Initialise the compare data
			s_compare_data.pp_char_in_seq[0] = &s_input_data.pp_char_in_seq[0][start_leaf_idx];
			s_compare_data.p_num_char_in_seq[0] = total_sequence_length;
			s_compare_data.pp_char_in_seq[1] = s_read_fasta.pp_char_in_seq[i];
			s_compare_data.p_num_char_in_seq[1] = s_read_fasta.p_num_char_in_seq[i];

			// printf("Sequence Length %d, %d\n", total_sequence_length, s_read_fasta.p_num_char_in_seq[i]);

			// Initialise the param config
			param_config_t s_param_config;
			s_param_config.match = 1; s_param_config.mismatch = -2;
			s_param_config.h = -5; s_param_config.g = -1;

			// Initialise the param config data
			output_stat_t s_output_stat;
			dp_state_t s_dp_state;
			int is_global = 0;
			find_alignment(&s_compare_data, &s_dp_state, &s_param_config, is_global, &s_output_stat);

			// if(s_output_stat.num_match < min_number_of_matches){
			// 	printf("Size of the read sequence %d, %d, %d\n", i, p_leaf_map_id[j], s_read_fasta.p_num_char_in_seq[i]);
			// 	// printf("Match less than exact match \n");
			// 	print_output_stat(&s_output_stat);

			// }
			// Calculating statistics:
			int align_length = s_output_stat.num_match + s_output_stat.num_mismatch + s_output_stat.num_gap;
			float percent_identity = 1.0*s_output_stat.num_match/align_length;
			float length_coverage = 1.0*align_length/s_read_fasta.p_num_char_in_seq[i];
			if(percent_identity >= 0.9 && length_coverage >= 0.8 && length_coverage > max_length_coverage){
				best_offset = s_compare_data.p_num_char_in_seq[0]+start_leaf_idx;
				max_length_coverage = length_coverage;
				max_percent_identity = percent_identity;
				best_align_length = align_length;
				// printf("Best Alignment offset %d, %d, %d, %d, %f, %f\n", i, j, best_offset, align_length, percent_identity, length_coverage);
			}
			// if(percent_identity >= max_percent_identity && length_coverage >= max_length_coverage){
			// 	best_offset = s_compare_data.p_num_char_in_seq[0]+start_leaf_idx;
			// 	best_align_length = align_length
			// 	max_length_coverage = length_coverage;
			// 	max_percent_identity = percent_identity;
			// 	// printf("Best Alignment offset %d, %d, %d, %d, %f, %f\n", i, j, best_offset, align_length, percent_identity, length_coverage);
			// }

			// printf("Best Alignment offset %d, %d, %d, %d, %d, %f, %f\n", i, j, s_compare_data.p_num_char_in_seq[0], start_leaf_idx, align_length, percent_identity, length_coverage);

			// printf("Offset alignment %d, %d\n", s_compare_data.p_num_char_in_seq[0], s_compare_data.p_num_char_in_seq[1]);
			// printf("Algnement value %d\n", s_compare_data.p_num_char_in_seq[0]+start_leaf_idx);

			// There are some cases where the Match value is less than 25.
			// Find out why this is happening

			// free memory for dp_state after compare

		}

		// if(max_length_coverage < 0.9){
		// 	printf("Error in max length coverage %d, %f, %f\n", i, max_length_coverage, max_percent_identity);
		// }
		if(max_length_coverage > 0.1){
			printf("READ_%d_from_Peach_reference.fasta %d %d\n", i+1, best_offset+1, best_offset+best_align_length+1);
		}
		else{
			printf("READ_%d_from_Peach_reference.fasta No hit found\n", i+1);
			num_misses += 1;
		}

		// free memory of s_read_fasta
		// for(int i=0;i<2;i++){
		// 	if(s_read_fasta.p_char_in_seq[i]){
		// 		free(s_read_fasta.p_char_in_seq[i]);
		// 	}
		// }
	}
	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

	printf("Num Search options %d, %d\n", num_local_search_options, num_misses);
	printf("Total time spent %f\n", time_spent);

	// free the allocated memory
	input_data_dealloc(&s_input_data, 1);
	input_data_dealloc(&s_compare_data, 0);
	input_data_dealloc(&s_read_fasta, num_read_seq);
	// for(int i=0;i<1;i++){
	// 	if(s_input_data.pp_char_in_seq[i]){
	// 		free(s_input_data.pp_char_in_seq[i]);
	// 	}
	// }
	// for(int i=0;i<num_read_seq;i++){
	// 	if(s_read_fasta.pp_char_in_seq[i]){
	// 		free(s_read_fasta.pp_char_in_seq[i]);
	// 	}
	// }
	// if(s_read_fasta.pp_char_in_seq){
	// 	free(s_read_fasta.pp_char_in_seq);
	// }
	// if(s_read_fasta.p_num_lines_in_seq){
	// 	free(s_read_fasta.p_num_lines_in_seq);
	// }
	// if(s_read_fasta.p_num_char_in_seq){
	// 	free(s_read_fasta.p_num_char_in_seq);
	// }
	if(p_leaf_map_id){
		free(p_leaf_map_id);
	}
	// free memory for the suffix tree
}