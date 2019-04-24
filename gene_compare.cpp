#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "fasta_read.h"
#include "gene_compare.h"
// Dumping output in required format
void print_output_stat(output_stat_t *p_output_stat){
		int total_seq_length = p_output_stat->num_match + \
			p_output_stat->num_mismatch + p_output_stat->num_gap;	
		float percent_gap = p_output_stat->num_gap*1.0/total_seq_length;
		float percent_identities = p_output_stat->num_match*1.0/total_seq_length;

		printf("Optimal Score %d\n", p_output_stat->opt_score);
		printf("Match %d, Mismatch %d, Gaps %d, Opening gap %d \n",
			p_output_stat->num_match, p_output_stat->num_mismatch, 
			p_output_stat->num_gap, p_output_stat->num_open_gap);
		printf("Precentage: Identities %f, Gaps %f\n", percent_identities,
			percent_gap);

		// // Dumping output sequence in aligned format with 60 characters per line
		// int seq_1_len = 0, seq_0_len = 0;
		// for(int i=total_seq_length; i>0;){
		// 	int width = 60;
		// 	if(width > i){
		// 		width = i;
		// 	}
		// 	printf("s1 %6d ", seq_1_len+1);
		// 	for(int j=0; j<width; j++){
		// 		printf("%c", p_output_stat->p_output_trace[2*(i-j-1)]);
		// 		if(p_output_stat->p_output_trace[2*(i-j-1)] != '-'){
		// 			seq_1_len++;
		// 		}
		// 	}
		// 	printf("%6d\n", seq_1_len+1);
		// 	printf("          ");
		// 	for(int j=0; j<width; j++){
		// 		if(p_output_stat->p_output_trace[2*(i-j-1)] ==
		// 		   p_output_stat->p_output_trace[2*(i-j-1)+1]){
		// 			printf("|");
		// 		}
		// 		else{
		// 			printf(" ");
		// 		}
		// 	}
		// 	printf("\n");
		// 	printf("s2 %6d ", seq_0_len+1);
		// 	for(int j=0; j<width; j++){
		// 		printf("%c", p_output_stat->p_output_trace[2*(i-j-1)+1]);
		// 		if(p_output_stat->p_output_trace[2*(i-j-1)+1] != '-'){
		// 			seq_0_len++;
		// 		}
		// 	}
		// 	printf("%6d\n", seq_0_len+1);
		// 	i -= width;
		// 	// printf("%d \n", i);
		// }
}

// Retrace the two sequences based on the input dynamic table
void retrace_global_old(input_data_t *p_input_data,
	dp_state_t *p_dp_state, int is_global,
	int local_idx_i, int local_idx_j,
	output_stat_t *p_output_stat,
	param_config_t *p_param_config){

	cell_value_t *p_dp_table = p_dp_state->p_dp_table;
	int seq_0_len = p_input_data->p_num_char_in_seq[0];
	int seq_1_len = p_input_data->p_num_char_in_seq[1];
	int i_idx = local_idx_i;
	int j_idx = local_idx_j;
	// int off = i_idx*(seq_0_len+1)+j_idx;

	// int s_cost = p_dp_table[off].cv[S_IDX];
	// int d_cost = p_dp_table[off].cv[D_IDX];
	// int i_cost = p_dp_table[off].cv[I_IDX];
	// // is_global == 0
	// max_t cur_max = find_max(s_cost, d_cost, i_cost, is_global);
	// int dir_idx = cur_max.idx;

	int num_match = 0, num_gap = 0, num_open_gap = 0, num_mismatch = 0;
	int is_prev_gap = 0, global_opt_score = -1*999999999;

	int next_dir_idx;
	int dir_idx;

	int cur_mat_pos = 0;

	int retrace_score = 0;
	// Initialise the cost
	int ms = p_param_config->match;
	int mms = p_param_config->mismatch;
	int hs = p_param_config->h;
	int gs = p_param_config->g;

	// while((i_idx > 0 && j_idx > 0) && dir_idx != -1){
	do{
		int off = i_idx*(seq_0_len+1)+j_idx;

		int s_cost = p_dp_table[off].cv[S_IDX];
		int d_cost = p_dp_table[off].cv[D_IDX];
		int i_cost = p_dp_table[off].cv[I_IDX];
		// is_global == 0
		max_t cur_max = find_max(s_cost, d_cost, i_cost, is_global);
		dir_idx = cur_max.idx;
		if(i_idx == local_idx_i && j_idx == local_idx_j){
			global_opt_score = cur_max.val;
			retrace_score = cur_max.val;
		}

		if(
			(retrace_score+hs != cur_max.val && is_prev_gap && dir_idx!=0) || 
		   (retrace_score != cur_max.val && (!is_prev_gap || (is_prev_gap && dir_idx==0))) ){
		   	int off_1 = (i_idx+1)*(seq_0_len+1)+j_idx+1;

			printf("Mismatching scores %d, %d, %d, %d, %d, %d\n", 
				i_idx, j_idx, retrace_score, cur_max.val, is_prev_gap, dir_idx);
			printf("Cost %d, %d\n", p_dp_table[off].cv[S_IDX], p_dp_table[off_1].cv[S_IDX]);
			printf("Cost %d, %d, %d\n", p_dp_table[off].cv[S_IDX], p_dp_table[off].cv[D_IDX],
				p_dp_table[off].cv[I_IDX]);

			printf("%c, %c, %c, %c, %c, %c\n", 
				p_input_data->pp_char_in_seq[1][i_idx-1],
				p_input_data->pp_char_in_seq[0][j_idx-1],
				p_input_data->pp_char_in_seq[1][i_idx],
				p_input_data->pp_char_in_seq[0][j_idx],
				p_input_data->pp_char_in_seq[1][i_idx+1],
				p_input_data->pp_char_in_seq[0][j_idx+1]);

			exit(0);
		}
		else{
			printf("Matching %d, %d, %d, %d\n", 
				retrace_score, cur_max.val, is_prev_gap, dir_idx);
			printf("Cost %d, %d, %d\n", p_dp_table[off].cv[S_IDX], p_dp_table[off].cv[D_IDX],
				p_dp_table[off].cv[I_IDX]);
			printf("Dir %d, %d, %d\n", p_dp_table[off].idx[S_IDX], p_dp_table[off].idx[D_IDX],
				p_dp_table[off].idx[I_IDX]);
		}

		next_dir_idx = p_dp_table[off].idx[dir_idx];		
		// printf("%d, %d, %d, %d, %d, %d, %d, %d, %d\n",
		// 	p_dp_table[off].cv[S_IDX],
		// 	p_dp_table[off].cv[D_IDX],
		// 	p_dp_table[off].cv[I_IDX],
		// 	p_dp_table[off].idx[S_IDX],
		// 	p_dp_table[off].idx[D_IDX],
		// 	p_dp_table[off].idx[I_IDX],
		// 	dir_idx, next_dir_idx, off);

		if(dir_idx == 0){ 
			// printf("%c%c %d,%d\n", \
			// p_input_data->p_char_in_seq[1][i_idx-1],\
			// p_input_data->p_char_in_seq[0][j_idx-1],
			// i_idx-1, j_idx-1);
			p_output_stat->p_output_trace[2*cur_mat_pos] = \
				p_input_data->pp_char_in_seq[1][i_idx-1];
			p_output_stat->p_output_trace[2*cur_mat_pos+1] = \
				p_input_data->pp_char_in_seq[0][j_idx-1];
			cur_mat_pos+=1;

			if(p_input_data->pp_char_in_seq[1][i_idx-1] == p_input_data->pp_char_in_seq[0][j_idx-1]){
				num_match += 1;
				retrace_score-=ms;
			}
			else{
				num_mismatch += 1;
				retrace_score-=mms;
			}
			is_prev_gap = 0;
		}
		// D
		else if(dir_idx == 1){ 
			// printf("%c- %d,%d\n", \
			// p_input_data->p_char_in_seq[1][i_idx-1],
			// i_idx-1, j_idx-1);
			p_output_stat->p_output_trace[2*cur_mat_pos] = \
				p_input_data->pp_char_in_seq[1][i_idx-1];
			p_output_stat->p_output_trace[2*cur_mat_pos+1] = '-';
			cur_mat_pos+=1;


			num_gap += 1;
			retrace_score-=gs;
			if(is_prev_gap == 0){
				num_open_gap += 1;
				retrace_score-=hs;
			}
			is_prev_gap = 1;
		}
		// I
		else if(dir_idx == 2){
			// printf("-%c %d,%d\n", \
			// p_input_data->p_char_in_seq[0][j_idx-1],
			// i_idx-1, j_idx-1);
			p_output_stat->p_output_trace[2*cur_mat_pos] = '-';
			p_output_stat->p_output_trace[2*cur_mat_pos+1] = \
				p_input_data->pp_char_in_seq[0][j_idx-1];
			cur_mat_pos+=1;

			num_gap += 1;
			retrace_score-=gs;
			if(is_prev_gap == 0){
				num_open_gap += 1;
				retrace_score-=hs;
			}
			is_prev_gap = 1;
		}
		// S 
		if(dir_idx == 0){ i_idx -= 1; j_idx -= 1;}
		// D
		else if(dir_idx == 1){ i_idx -= 1; }
		// I
		else if(dir_idx == 2){ j_idx -= 1; }


		// Update
		// off = i_idx*(seq_0_len+1)+j_idx;
		// dir_idx = next_dir_idx;
	}while((i_idx > 0 && j_idx > 0) && dir_idx != -1);

	p_output_stat->num_match = num_match;
	p_output_stat->num_mismatch = num_mismatch;
	p_output_stat->num_gap = num_gap;
	p_output_stat->num_open_gap = num_open_gap;
	p_output_stat->opt_score = global_opt_score;
	if(!is_global){
		p_input_data->p_num_char_in_seq[0] = j_idx-1;
		p_input_data->p_num_char_in_seq[1] = i_idx-1;
		// printf("i, j idx %d, %d\n", i_idx, j_idx);
	}
}






// Retrace the two sequences based on the input dynamic table
void retrace_global(input_data_t *p_input_data,
	dp_state_t *p_dp_state, int is_global,
	int local_idx_i, int local_idx_j,
	output_stat_t *p_output_stat,
	param_config_t *p_param_config){

	cell_value_t *p_dp_table = p_dp_state->p_dp_table;
	int seq_0_len = p_input_data->p_num_char_in_seq[0];
	int seq_1_len = p_input_data->p_num_char_in_seq[1];
	int i_idx = local_idx_i;
	int j_idx = local_idx_j;
	int off = i_idx*(seq_0_len+1)+j_idx;

	int s_cost = p_dp_table[off].cv[S_IDX];
	int d_cost = p_dp_table[off].cv[D_IDX];
	int i_cost = p_dp_table[off].cv[I_IDX];
	// is_global == 0
	max_t cur_max = find_max(s_cost, d_cost, i_cost, is_global);
	int dir_idx = cur_max.idx;


	int num_match = 0, num_gap = 0, num_open_gap = 0, num_mismatch = 0;
	int is_prev_gap = 0, global_opt_score = cur_max.val;//-1*999999999;

	int next_dir_idx;
	// int dir_idx;

	int cur_mat_pos = 0;

	int retrace_score = cur_max.val;
	// Initialise the cost
	int ms = p_param_config->match;
	int mms = p_param_config->mismatch;
	int hs = p_param_config->h;
	int gs = p_param_config->g;

	while((i_idx > 0 && j_idx > 0) && dir_idx != -1){
		// Filling the S,D,I values for each cell
		int i_off = i_idx*(seq_0_len+1), i_1_off = (i_idx-1)*(seq_0_len+1);
		int j_off = j_idx, j_1_off = (j_idx-1);
		int i_1_j_off = i_1_off+j_off, i_j_1_off = i_off + j_1_off;
		int i_1_j_1_off = i_1_off + j_1_off, i_j_off = i_off + j_off;

		// if(
		// 	(retrace_score+hs+gs != cur_max.val && is_prev_gap) || 
		//    (retrace_score != cur_max.val && !is_prev_gap) ){
		//    	int off_1 = (i_idx+1)*(seq_0_len+1)+j_idx+1;

		// 	printf("Mismatching scores %d, %d, %d, %d, %d, %d\n", 
		// 		i_idx, j_idx, retrace_score, cur_max.val, is_prev_gap, dir_idx);
		// 	// printf("Cost %d, %d\n", p_dp_table[off].cv[S_IDX], p_dp_table[off_1].cv[S_IDX]);
		// 	printf("Cost %d, %d, %d\n", p_dp_table[i_j_off].cv[S_IDX], p_dp_table[i_j_off].cv[D_IDX],
		// 		p_dp_table[i_j_off].cv[I_IDX]);

		// 	// printf("%c, %c, %c, %c, %c, %c\n", 
		// 	// 	p_input_data->p_char_in_seq[1][i_idx-1],
		// 	// 	p_input_data->p_char_in_seq[0][j_idx-1],
		// 	// 	p_input_data->p_char_in_seq[1][i_idx],
		// 	// 	p_input_data->p_char_in_seq[0][j_idx],
		// 	// 	p_input_data->p_char_in_seq[1][i_idx+1],
		// 	// 	p_input_data->p_char_in_seq[0][j_idx+1]);

		// 	exit(0);
		// }
		// else{
		// 	printf("Matching %d, %d, %d, %d\n", 
		// 		retrace_score, cur_max.val, is_prev_gap, dir_idx);
		// 	printf("Cost %d, %d, %d\n", p_dp_table[i_j_off].cv[S_IDX], p_dp_table[i_j_off].cv[D_IDX],
		// 		p_dp_table[i_j_off].cv[I_IDX]);
		// 	printf("Dir %d, %d, %d\n", p_dp_table[i_j_off].idx[S_IDX], p_dp_table[i_j_off].idx[D_IDX],
		// 		p_dp_table[i_j_off].idx[I_IDX]);
		// }

		if(dir_idx == 0){ 
			// printf("%c%c %d,%d\n", \
			// p_input_data->p_char_in_seq[1][i_idx-1],\
			// p_input_data->p_char_in_seq[0][j_idx-1],
			// i_idx-1, j_idx-1);
			p_output_stat->p_output_trace[2*cur_mat_pos] = \
				p_input_data->pp_char_in_seq[1][i_idx-1];
			p_output_stat->p_output_trace[2*cur_mat_pos+1] = \
				p_input_data->pp_char_in_seq[0][j_idx-1];
			cur_mat_pos+=1;

			if(p_input_data->pp_char_in_seq[1][i_idx-1] == p_input_data->pp_char_in_seq[0][j_idx-1]){
				num_match += 1;
				retrace_score-=ms;
			}
			else{
				num_mismatch += 1;
				retrace_score-=mms;
			}
			is_prev_gap = 0;

			i_idx -= 1; j_idx -= 1;


			// Find the max from the (i-1,j-1) cell
			cur_max = find_max(p_dp_table[i_1_j_1_off].cv[S_IDX],
				p_dp_table[i_1_j_1_off].cv[D_IDX],
				p_dp_table[i_1_j_1_off].cv[I_IDX], is_global);
			dir_idx = cur_max.idx;
		}
		// D
		else if(dir_idx == 1){ 
			// printf("%c- %d,%d\n", \
			// p_input_data->p_char_in_seq[1][i_idx-1],
			// i_idx-1, j_idx-1);
			p_output_stat->p_output_trace[2*cur_mat_pos] = \
				p_input_data->pp_char_in_seq[1][i_idx-1];
			p_output_stat->p_output_trace[2*cur_mat_pos+1] = '-';
			cur_mat_pos+=1;


			num_gap += 1;
			retrace_score-=gs;
			if(is_prev_gap == 0){
				num_open_gap += 1;
				retrace_score-=hs;
			}
			is_prev_gap = 1;

			i_idx -= 1; 


			// Filling the D value based on the previous values
			int s_cost = p_dp_table[i_1_j_off].cv[S_IDX] + hs + gs;
			int d_cost = p_dp_table[i_1_j_off].cv[D_IDX] + gs;
			int i_cost = p_dp_table[i_1_j_off].cv[I_IDX] + hs + gs;

			cur_max = find_max(s_cost, d_cost, i_cost, is_global);
			dir_idx = cur_max.idx;
		}
		// I
		else if(dir_idx == 2){
			// printf("-%c %d,%d\n", \
			// p_input_data->p_char_in_seq[0][j_idx-1],
			// i_idx-1, j_idx-1);
			p_output_stat->p_output_trace[2*cur_mat_pos] = '-';
			p_output_stat->p_output_trace[2*cur_mat_pos+1] = \
				p_input_data->pp_char_in_seq[0][j_idx-1];
			cur_mat_pos+=1;

			num_gap += 1;
			retrace_score-=gs;
			if(is_prev_gap == 0){
				num_open_gap += 1;
				retrace_score-=hs;
			}
			is_prev_gap = 1;

			j_idx -= 1;


			// Filling the I value based on the previous values
			s_cost = p_dp_table[i_j_1_off].cv[S_IDX] + hs + gs;
			d_cost = p_dp_table[i_j_1_off].cv[D_IDX] + hs + gs;
			i_cost = p_dp_table[i_j_1_off].cv[I_IDX] + gs;

			cur_max = find_max(s_cost, d_cost, i_cost, is_global);
			dir_idx = cur_max.idx;
		}
	}

	p_output_stat->num_match = num_match;
	p_output_stat->num_mismatch = num_mismatch;
	p_output_stat->num_gap = num_gap;
	p_output_stat->num_open_gap = num_open_gap;
	p_output_stat->opt_score = global_opt_score;
	if(!is_global){
		p_input_data->p_num_char_in_seq[0] = j_idx-1;
		p_input_data->p_num_char_in_seq[1] = i_idx-1;
		// printf("i, j idx %d, %d\n", i_idx, j_idx);
	}
}



// Function to align the two input sequences
void find_alignment(input_data_t *p_input_data, 
	dp_state_t *p_dp_state,
	param_config_t *p_param_config,
	int is_global,
	output_stat_t *p_output_stat){

	int seq_0_len = p_input_data->p_num_char_in_seq[0];
	int seq_1_len = p_input_data->p_num_char_in_seq[1];
	int table_size = (seq_0_len+1)*(seq_1_len+1);
	int max_length = seq_0_len+seq_1_len;

	// Based on the input dimensions allocate the necessary memory required for alignment
	p_dp_state->p_dp_table = (cell_value_t *)malloc(table_size*sizeof(cell_value_t));
	if(!p_dp_state->p_dp_table){
		printf("Error in allocatin memory for the table: %d\n", table_size);
		exit(-1);
	}
	p_output_stat->p_output_trace = (char *)malloc(3*max_length*sizeof(char));
	if(!p_output_stat->p_output_trace){
		printf("Error in allocatin memory for output trace: %d\n", 3*max_length);
		exit(-1);
	}
	cell_value_t *p_dp_table = p_dp_state->p_dp_table;
	char *p_seq_0 = p_input_data->pp_char_in_seq[0];
	char *p_seq_1 = p_input_data->pp_char_in_seq[1];

	// Initialise the cost
	int ms = p_param_config->match;
	int mms = p_param_config->mismatch;
	int hs = p_param_config->h;
	int gs = p_param_config->g;

	int max_local_cost = -1*99999999;
	int local_idx_i, local_idx_j;

	// printf("Cost %d, %d, %d, %d \n", ms, mms, hs, gs);

	// Initialise the init values with fixed numbers
	for(int i=0; i<3; i++) {
		p_dp_table[0].cv[i] = 0;
		p_dp_table[0].idx[i] = 0;
	}
	for (int i=1; i<=seq_0_len; i++){
		p_dp_table[i].cv[S_IDX] = hs+gs*seq_0_len-1;
		p_dp_table[i].cv[D_IDX] = hs+gs*seq_0_len-1;
		p_dp_table[i].cv[I_IDX] = hs+gs*i;
		for(int j=0; j<3; j++) p_dp_table[i].idx[j] = 2;
	}
	for (int i=1; i<=seq_1_len; i++){
		p_dp_table[i*(seq_0_len+1)].cv[S_IDX] = hs+gs*seq_1_len-1;
		p_dp_table[i*(seq_0_len+1)].cv[I_IDX] = hs+gs*seq_1_len-1;
		p_dp_table[i*(seq_0_len+1)].cv[D_IDX] = hs+gs*i;
		for(int j=0; j<3; j++) p_dp_table[i*(seq_0_len+1)].idx[j] = 1;
	}
	// Calculation of the full table
	for (int i=1; i<=seq_1_len; i++){
		for(int j=1; j<=seq_0_len; j++){
			// Filling the S,D,I values for each cell
			int i_off = i*(seq_0_len+1), i_1_off = (i-1)*(seq_0_len+1);
			int j_off = j, j_1_off = (j-1);
			int i_1_j_off = i_1_off+j_off, i_j_1_off = i_off + j_1_off;
			int i_1_j_1_off = i_1_off + j_1_off, i_j_off = i_off + j_off;
			// Filling the S value based on the previous values
			int cost = mms;
			// Check the match for the current cell
			if(p_seq_0[j-1] == p_seq_1[i-1]) cost = ms;
			// Find the max from the (i-1,j-1) cell
			max_t prev_max = find_max(p_dp_table[i_1_j_1_off].cv[S_IDX],
				p_dp_table[i_1_j_1_off].cv[D_IDX],
				p_dp_table[i_1_j_1_off].cv[I_IDX], is_global);
			p_dp_table[i_j_off].cv[S_IDX] = prev_max.val + cost;
			p_dp_table[i_j_off].idx[S_IDX] = prev_max.idx;


			if(p_dp_table[i_j_off].cv[S_IDX] >= max_local_cost){
				max_local_cost = p_dp_table[i_j_off].cv[S_IDX];
				local_idx_i = i; local_idx_j = j;
				// printf("The cost %d, %d, %d, %d\n", i, j, i_j_off, p_dp_table[i_j_off].cv[S_IDX]);
			}

			// Filling the D value based on the previous values
			int s_cost = p_dp_table[i_1_j_off].cv[S_IDX] + hs + gs;
			int d_cost = p_dp_table[i_1_j_off].cv[D_IDX] + gs;
			int i_cost = p_dp_table[i_1_j_off].cv[I_IDX] + hs + gs;

			// if(i==9675+1 && j == 9841){
			// 	printf("Forward Pass %d, %d\n", i, j);
			// 	printf("%d, %d, %d\n", s_cost, d_cost, i_cost);
			// 	printf("%d, %d, %d\n", p_dp_table[i_1_j_off].cv[S_IDX], 
			// 		p_dp_table[i_1_j_off].cv[D_IDX],
			// 		p_dp_table[i_1_j_off].cv[I_IDX]);
			// }

			max_t cur_max = find_max(s_cost, d_cost, i_cost, is_global);
			p_dp_table[i_j_off].cv[D_IDX] = cur_max.val;
			p_dp_table[i_j_off].idx[D_IDX] = cur_max.idx;

			// Filling the I value based on the previous values
			s_cost = p_dp_table[i_j_1_off].cv[S_IDX] + hs + gs;
			d_cost = p_dp_table[i_j_1_off].cv[D_IDX] + hs + gs;
			i_cost = p_dp_table[i_j_1_off].cv[I_IDX] + gs;

			cur_max = find_max(s_cost, d_cost, i_cost, is_global);
			p_dp_table[i_j_off].cv[I_IDX] = cur_max.val;
			p_dp_table[i_j_off].idx[I_IDX] = cur_max.idx;

			// printf("%d, %d: ", i, j);
			// for(int k=0; k<3; k++) printf("%d ", p_dp_table[i_j_off].cv[k]);
			// // printf("\n");
			// for(int k=0; k<3; k++) printf("%d ", p_dp_table[i_j_off].idx[k]);
			// printf("\n");
		}
	}
	// Let the retrace begin
	// output_stat_t output_stat;
	if(is_global) {
		retrace_global(p_input_data, p_dp_state, is_global,
			seq_1_len, seq_0_len, p_output_stat, p_param_config);
	}
	else {
		retrace_global(p_input_data, p_dp_state, is_global,
			local_idx_i, local_idx_j, p_output_stat, p_param_config);
	}

	if(p_dp_state->p_dp_table){
		free(p_dp_state->p_dp_table);
	}
	if(p_output_stat->p_output_trace){
		free(p_output_stat->p_output_trace);
	}
}
