#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>

#include "fasta_read.h"
#include "suffix_tree.h"

// #define __PRINTF__(x) printf(x);
// #define __PRINTF__(x)

void __PRINTF__(const char *format, ...)
{

 //    va_list args;
 //    va_start(args, format);
	// vprintf(format, args);
 //    va_end(args);
}

// typedef struct _st_node_t_ {
// 	int node_idx; // Unique id for each node
// 	int leaf_idx; // Unique id for each suffix. All leaf less node had index -1
// 	struct _st_node_t_ *p_suffix;
// 	struct _st_node_t_ *p_child;
// 	int edge_label[2]; // start and end offset of the edge
// 	int string_depth; // Offset to the string from the root
// 	struct _st_node_t_ *p_sibling;
// }st_node_t;

// void find_node(st_node_t *p_node, char *p_char_in_seq, int offset){
// 	// 
// 	int start_idx = p_node->edge_label[1] - p_node->string_depth;
// 	int current_string_offset = offset + p_node->string_depth;
// 	st_node_t *p_child = p_node->p_child;

// 	if(p_child == 0){
// 		return p_node;
// 	}
// 	// Search all the siblings until one matches a label
// 	while(p_child){
// 		if(p_char_in_seq[])
// 		// go to the next child
// 		p_child = p_child->p_sibling;
// 	}
// }

// Logic used:
// https://eecs.wsu.edu/~ananth/CptS571/Lectures/iScribe_STconstruction_McCreightsAlgo.pdf

int GLOBAL_NODE_IDX = 0;
// Wrapper function for initialising tha node
st_node_t *create_node(int node_idx, int leaf_idx,
	st_node_t *p_suffix, st_node_t *p_child, int lbl_start_idx,
	int lbl_end_idx, int parent_string_depth, st_node_t *p_sibling){

	// Add the full string as the first suffix
	st_node_t *p_node = (st_node_t *)malloc(1*sizeof(st_node_t));
	if(p_node == 0){
		printf("Error in allocating input node\n");
		exit(-1);
	}
	// Initiliase node
	p_node->node_idx = GLOBAL_NODE_IDX;//node_idx;
	GLOBAL_NODE_IDX += 1;
	p_node->leaf_idx = leaf_idx;
	p_node->p_suffix = p_suffix;
	p_node->p_child = p_child;
	p_node->edge_label[0] = lbl_start_idx;
	p_node->edge_label[1] = lbl_end_idx;
	p_node->string_depth = parent_string_depth + lbl_end_idx - lbl_start_idx;
	p_node->p_sibling = p_sibling;
	p_node->leaf_idx_range[0] = -1;
	p_node->leaf_idx_range[1] = -1;
	return p_node;
}

void fill_suffix_link_info(suffix_link_info_t *p_sl_info,
	st_node_t *p_u_node, st_node_t *p_u_dash_node){

	p_sl_info->p_u_dash_node = p_u_dash_node;
	p_sl_info->p_u_node = p_u_node;
	if(p_u_dash_node){
		__PRINTF__("SL info update %d, %d\n", p_u_node->node_idx, p_u_dash_node->node_idx);
	}
	else{
		__PRINTF__("SL info update %d\n", p_u_node->node_idx);
	}
	// p_sl_info->beta = beta;
}

void init_leaf_idx_range(st_node_t *p_root_node, int *leaf_idx_map, int *cur_arr_idx){
	st_node_t *p_parent_node = p_root_node;
	if(p_parent_node){
		// If the current node is a leaf node initialise the leaf index mapping
		if(p_parent_node->leaf_idx != -1){

			leaf_idx_map[cur_arr_idx[0]] = p_parent_node->leaf_idx;

			if(1){//p_parent_node->string_depth >= 25){
				p_parent_node->leaf_idx_range[0] = cur_arr_idx[0];
				p_parent_node->leaf_idx_range[1] = cur_arr_idx[0];
			}
			
			// printf("Inc %d, ", cur_arr_idx[0]);
			cur_arr_idx[0]++;
		}
		// Before calling the child initialise the left index [0]
		// printf("0 %d, ", cur_arr_idx[0]);

		init_leaf_idx_range(p_parent_node->p_child, leaf_idx_map, cur_arr_idx);
		init_leaf_idx_range(p_parent_node->p_sibling, leaf_idx_map, cur_arr_idx);

		// After processing all the children initialise left and right nodes
		if(p_parent_node->leaf_idx == -1){
			p_parent_node->leaf_idx_range[0] = p_parent_node->p_child->leaf_idx_range[0];
			st_node_t *p_sib_node = p_parent_node->p_child;
			// printf("child %d\n", p_parent_node->p_child);
			while(p_sib_node->p_sibling){
				p_sib_node = p_sib_node->p_sibling;
			}
			// printf("sibling %d\n", p_sib_node);
			p_parent_node->leaf_idx_range[1] = p_sib_node->leaf_idx_range[1];
		}

	}
}

void dfs(st_node_t *p_root_node, input_data_t *p_input_data, int *p_max_depth){
	char *p_char_in_seq = p_input_data->pp_char_in_seq[0];
	int num_char_in_seq = p_input_data->p_num_char_in_seq[0];

	st_node_t *p_parent_node = p_root_node;
	if(p_parent_node){
		dfs(p_parent_node->p_child, p_input_data, p_max_depth);
		if(p_parent_node->leaf_idx != -1){
			// printf("Node idx %d, %d, \n",p_parent_node->node_idx,
			// 	p_parent_node->leaf_idx);
			int bwt_idx = p_parent_node->leaf_idx-1;
			if(bwt_idx<0){
				bwt_idx = num_char_in_seq-1;
			}
			// printf("%c, %d\n",p_char_in_seq[bwt_idx], p_parent_node->leaf_idx);
			//printf("%c\n",p_char_in_seq[bwt_idx]);
		}

		if(p_parent_node->leaf_idx == -1){
			if(p_max_depth[0] < p_parent_node->string_depth){
				p_max_depth[0] = p_parent_node->string_depth;
			}
		}

		// if(p_parent_node->node_idx == 1883 || p_parent_node->node_idx == 2114
		// || p_parent_node->node_idx == 1914 || p_parent_node->node_idx == 179
		// || p_parent_node->node_idx == 93 || p_parent_node->node_idx == 181
		// || p_parent_node->node_idx == 111
		// || p_parent_node->node_idx == 156 || p_parent_node->node_idx == 157){

		// if(p_parent_node->node_idx == 104 || p_parent_node->node_idx == 46
		// || p_parent_node->node_idx == 17  || p_parent_node->node_idx == 21
		// || p_parent_node->node_idx == 99  || p_parent_node->node_idx == 4){

		if(0){

			int suf_node_idx = -1;
			if(p_parent_node->p_suffix){
				suf_node_idx = p_parent_node->p_suffix->node_idx;
			}
			printf("Node idx %d, %d, %d, %d [",p_parent_node->node_idx,
					p_parent_node->leaf_idx, p_parent_node->string_depth,
					suf_node_idx);

			int end_idx = p_parent_node->edge_label[1];
			if(p_parent_node->edge_label[1] > p_parent_node->edge_label[0] + 10){
				end_idx = p_parent_node->edge_label[0] + 10;
			}
			for(int i=p_parent_node->edge_label[0] ; i < end_idx; i++){
				printf("%c ", p_char_in_seq[i]);
			}
			printf("] [");
			st_node_t *p_child = p_parent_node->p_child;
			while(p_child){
				printf("%d ",p_child->node_idx);
				p_child = p_child->p_sibling;
			}
			printf("]\n");
		}

		dfs(p_parent_node->p_sibling, p_input_data, p_max_depth);
	}
}

st_node_t * find_node_matching_path(st_node_t *p_root_node, char *ref_seq,
	char *p_char_in_seq, int leaf_offset, int size_of_seq, int min_number_of_matches){
	// st_node_t *p_parent_node = p_root_node;
	// st_node_t *p_prev_parent_node = p_sl_info->p_v_dash_node;
	// printf("Find path %d %d %d\n", leaf_offset, size_of_seq, min_number_of_matches);

	if(size_of_seq-min_number_of_matches < 0){
		printf("This condition cannot happen %d, %d\n", size_of_seq, min_number_of_matches);
		exit(-1);
	}
	st_node_t *max_match_depth_node = 0;
	int max_num_matches = 0;
	for(int j=leaf_offset; j<=size_of_seq-min_number_of_matches; j++){
		int num_matches = 0;
		st_node_t *p_parent_node = p_root_node;
		// st_node_t *p_prev_parent_node = p_root_node;
		int cur_char_offset = j;
		// printf("Matching offset %d, %d\n", j, size_of_seq-min_number_of_matches);

		// Iterate the depth of the tree from root to leaf
		while(p_parent_node){
			// Compare the edge label of each children starting from the leftmost sibling
			st_node_t *p_sib_node = p_parent_node->p_child;
			st_node_t *p_prev_sib_node = 0;
			// Iterate through all the sibling till you find a match
			int is_match_found = 0;
			int partial_match_offset = 0, start_idx, end_idx;
			while(p_sib_node){
				start_idx = p_sib_node->edge_label[0];
				end_idx = p_sib_node->edge_label[1];
				// printf("Comparing values %c, %c\n", p_char_in_seq[cur_char_offset], ref_seq[start_idx]);
				// If one character matches then node exist
				if(p_char_in_seq[cur_char_offset] == ref_seq[start_idx]){
					is_match_found = 1;
					for(int i = 0; i < end_idx-start_idx ; i++){
						if(p_char_in_seq[cur_char_offset+i] != ref_seq[start_idx+i]){
							// It is a partial match
							partial_match_offset = i;
							is_match_found = 2;
							break;
						}
					}
					break;
				}
				// if the character did not match and character is lesser than current value
				// it is a mismatch
				else if(p_char_in_seq[cur_char_offset] <= ref_seq[start_idx]){
					break;
				}
				p_prev_sib_node = p_sib_node;
				p_sib_node = p_sib_node->p_sibling;
			}
			// is_match_found = 1 : Its a perfect match and so need to goto the next node
			if(is_match_found == 1){
				// Make the current sibling as the next parent and go down the tree
				// p_prev_parent_node = p_parent_node;
				p_parent_node = p_sib_node;
				cur_char_offset += (end_idx-start_idx);
				num_matches += (end_idx-start_idx);
				if(num_matches != p_parent_node->string_depth){
					printf("Error in string depth %d, %d\n", num_matches, p_parent_node->string_depth);
				}
				// printf("Error in string depth %d, %d\n", num_matches, p_parent_node->string_depth);

				// printf("matched value %d, %d, %d\n", j, num_matches, p_parent_node->node_idx);
			}
			else{
				num_matches += partial_match_offset;
				// st_node_t *p_matching_node;
				// if(is_match_found == 0) p_matching_node = p_parent_node;
				// else p_matching_node = p_prev_parent_node;

				// printf("partial matched value %d, %d\n", j, num_matches);
				if(num_matches > max_num_matches && num_matches >= min_number_of_matches){
				// if(num_matches > max_num_matches && p_parent_node->string_depth >= min_number_of_matches){

					// if(is_match_found == 0) max_match_depth_node = p_parent_node;
					// else max_match_depth_node = p_prev_parent_node;
					max_match_depth_node = p_parent_node;

					// printf("The matching %d, %d, %d, %d, %d, %d\n", max_match_depth_node->leaf_idx_range[0],
					// 	max_match_depth_node->leaf_idx_range[1], num_matches, is_match_found,
					// 	max_match_depth_node->string_depth, partial_match_offset);

					max_num_matches = num_matches;
				}
				break;
			}
		}
	}

	// printf("Final return %d, %d\n", max_num_matches, max_match_depth_node);
	return max_match_depth_node;
}

st_node_t * find_node_matching_path_sl(st_node_t *p_root_node, char *ref_seq,
	char *p_char_in_seq, int leaf_offset, int size_of_seq, int min_number_of_matches){
	// st_node_t *p_parent_node = p_root_node;
	// st_node_t *p_prev_parent_node = p_sl_info->p_v_dash_node;
	// printf("Find path %d %d %d\n", leaf_offset, size_of_seq, min_number_of_matches);

	if(size_of_seq-min_number_of_matches < 0){
		printf("This condition cannot happen %d, %d\n", size_of_seq, min_number_of_matches);
		exit(-1);
	}
	st_node_t *max_match_depth_node = 0;
	int max_num_matches = 0;
	int num_matches = 0;
	st_node_t *p_parent_node = p_root_node;
	st_node_t *p_prev_parent_node = p_root_node;
	int cur_char_offset = leaf_offset;
	for(int j=leaf_offset; j<=size_of_seq-min_number_of_matches; j++){
		
		// printf("Matching offset %d, %d\n", j, size_of_seq-min_number_of_matches);

		// Iterate the depth of the tree from root to leaf
		while(p_parent_node){
			// Compare the edge label of each children starting from the leftmost sibling
			st_node_t *p_sib_node = p_parent_node->p_child;
			st_node_t *p_prev_sib_node = 0;
			// Iterate through all the sibling till you find a match
			int is_match_found = 0;
			int partial_match_offset = 0, start_idx, end_idx;
			while(p_sib_node){
				start_idx = p_sib_node->edge_label[0];
				end_idx = p_sib_node->edge_label[1];
				// printf("Comparing values %c, %c\n", p_char_in_seq[cur_char_offset], ref_seq[start_idx]);
				// If one character matches then node exist
				if(p_char_in_seq[cur_char_offset] == ref_seq[start_idx]){
					is_match_found = 1;
					for(int i = 0; i < end_idx-start_idx ; i++){
						if(p_char_in_seq[cur_char_offset+i] != ref_seq[start_idx+i]){
							// It is a partial match
							partial_match_offset = i;
							is_match_found = 2;
							break;
						}
					}
					break;
				}
				// if the character did not match and character is lesser than current value
				// it is a mismatch
				else if(p_char_in_seq[cur_char_offset] <= ref_seq[start_idx]){
					break;
				}
				p_prev_sib_node = p_sib_node;
				p_sib_node = p_sib_node->p_sibling;
			}
			// is_match_found = 1 : Its a perfect match and so need to goto the next node
			if(is_match_found == 1){
				// Make the current sibling as the next parent and go down the tree
				p_prev_parent_node = p_parent_node;
				p_parent_node = p_sib_node;
				cur_char_offset += (end_idx-start_idx);
				num_matches += (end_idx-start_idx);
				// printf("matched value %d, %d, %d\n", j, num_matches, p_parent_node->node_idx);
			}
			else{
				num_matches += partial_match_offset;
				// printf("partial matched value %d, %d, %d\n", j, num_matches, p_parent_node->node_idx);
				if(num_matches > max_num_matches && num_matches >= min_number_of_matches){

					if(is_match_found == 0) max_match_depth_node = p_parent_node;
					else max_match_depth_node = p_prev_parent_node;

					// printf("The matching %d, %d, %d, %d\n", p_parent_node->leaf_idx_range[0],
					// 	p_parent_node->leaf_idx_range[1], num_matches, is_match_found);

					max_num_matches = num_matches;
				}

				// Initialise the values for the next offset
				num_matches -= (partial_match_offset+1);
				if(is_match_found == 0) {
					p_parent_node = p_parent_node->p_suffix;
					p_prev_parent_node = p_parent_node;
				}
				else {
					p_parent_node = p_prev_parent_node->p_suffix;
					p_prev_parent_node = p_parent_node;
				}
				// cur_char_offset += 1;
				if(p_parent_node->node_idx == 0){
					cur_char_offset = j+1;
					num_matches = 0;
				}
				// printf("next offset %d, %d, %d\n", cur_char_offset, num_matches, p_parent_node->node_idx);
				break;
			}
		}
	}

	// printf("Final return %d, %d\n", max_num_matches, max_match_depth_node);
	return max_match_depth_node;
}


void find_path(st_node_t *p_root_node, char *p_char_in_seq, int leaf_offset, int size_of_seq,
	suffix_link_info_t *p_sl_info, int leaf_idx){

	st_node_t *p_parent_node = p_root_node;
	st_node_t *p_prev_parent_node = p_sl_info->p_v_dash_node;
	int cur_char_offset = leaf_offset;
	__PRINTF__("Find path %d %d\n", p_root_node->node_idx, leaf_offset);

	// Iterate the depth of the tree from root to leaf
	while(p_parent_node){
		// Compare the edge label of each children starting from the leftmost sibling
		st_node_t *p_sib_node = p_parent_node->p_child;
		st_node_t *p_prev_sib_node = 0;
		// Iterate through all the sibling till you find a match
		int is_match_found = 0;
		int partial_match_offset = -1, start_idx, end_idx;
		while(p_sib_node){
			start_idx = p_sib_node->edge_label[0];
			end_idx = p_sib_node->edge_label[1];
			// If one character matches then node exist
			if(p_char_in_seq[cur_char_offset] == p_char_in_seq[start_idx]){
				is_match_found = 1;
				for(int i = 0; i < end_idx-start_idx ; i++){
					if(p_char_in_seq[cur_char_offset+i] != p_char_in_seq[start_idx+i]){
						// It is a partial match
						partial_match_offset = i;
						is_match_found = 2;
						break;
					}
				}
				break;
			}
			// if the character did not match and character is lesser than current value
			// it is a mismatch
			else if(p_char_in_seq[cur_char_offset] <= p_char_in_seq[start_idx]){
				break;
			}
			p_prev_sib_node = p_sib_node;
			p_sib_node = p_sib_node->p_sibling;
		}
		// Based on whether the siblings match or not we have three cases
		// is_match_found = 0 : If no match is found, add this character sequence as a new node
		if(is_match_found == 0){
			// Add the character sequence into current parent node
			// The new node created is a leaf node and leaf node uses the offset as the index
			st_node_t *p_new_node = create_node(-1, leaf_idx, (st_node_t *)0 /*TBD*/, 
				(st_node_t *)0 /*child */, 
				cur_char_offset, size_of_seq, p_parent_node->string_depth, 
				p_sib_node);
			if(p_prev_sib_node){
				// if p_prev_sib_node is not null then add the current node as the next sibling
				p_prev_sib_node->p_sibling = p_new_node;
			}
			else{
				// if p_prev_sib_node is null then add the current node as the leftmost node
				// Make the new node the first sibling
				p_parent_node->p_child = p_new_node;
				// This condition doesnot happen other than the root node case since other cases would always split and generate two node minimum
			}
			// printf("Create node %d\n", p_new_node->node_idx);

			// All cases should result in the parent knowing suffix links
			// Since either it would be aroot node or internal which would have been created earlier
			if(p_parent_node->p_suffix){
				fill_suffix_link_info(p_sl_info, p_parent_node, (st_node_t *)0);
			}
			else{
				fill_suffix_link_info(p_sl_info, p_parent_node, p_prev_parent_node);
				// printf("This Condition should not happen %d\n", p_parent_node->node_idx);
				// exit(-1);
			}
			break;
		}
		// is_match_found = 2 : Its a partial match. So break and introduce a new node
		if(is_match_found == 2){
			// printf("Current node matching %d %d 0x%x %d\n", is_match_found, partial_match_offset,
			// 	p_sib_node,
			// 	cur_char_offset + partial_match_offset);
			int parent_string_depth = p_sib_node->string_depth - (p_sib_node->edge_label[1] - p_sib_node->edge_label[0] - partial_match_offset);
			st_node_t *p_new_node = create_node(-1, -1, (st_node_t *)0, 
				p_sib_node /*child */, 
				p_sib_node->edge_label[0], 
				p_sib_node->edge_label[0] + partial_match_offset, 
				p_parent_node->string_depth, 
				p_sib_node->p_sibling/* sibling */);

			// Update the sibling with the partial node information
			p_sib_node->edge_label[0] += partial_match_offset;
			p_sib_node->p_sibling = (st_node_t *)0;

			st_node_t *p_new_leaf_node = create_node(-1, leaf_idx, (st_node_t *)0 /*TBD*/, 
				(st_node_t *)0 /*child */, 
				cur_char_offset + partial_match_offset, size_of_seq, 
				parent_string_depth, 
				(st_node_t *)0 /*Update this after comparing between the two new nodes */);

			if(p_prev_sib_node){
				p_prev_sib_node->p_sibling = p_new_node;
			}
			else{
				p_parent_node->p_child = p_new_node; 
			}

			if(p_char_in_seq[p_sib_node->edge_label[0]] < p_char_in_seq[p_new_leaf_node->edge_label[0]]){
			   	p_new_node->p_child = p_sib_node;
				p_sib_node->p_sibling = p_new_leaf_node;
				p_new_leaf_node->p_sibling = (st_node_t *)0;
			}
			else{
				p_new_node->p_child = p_new_leaf_node;
				p_new_leaf_node->p_sibling = p_sib_node;
				p_sib_node->p_sibling = (st_node_t *)0;
			}
			// Update the suffix link info
			fill_suffix_link_info(p_sl_info, p_new_node, p_parent_node);
			break;
		}

		// is_match_found = 1 : Its a perfect match and so need to goto the next node
		if(is_match_found == 1){
			// Make the current sibling as the next parent and go down the tree
			p_prev_parent_node = p_parent_node;
			p_parent_node = p_sib_node;
			cur_char_offset += (end_idx-start_idx);
		}
	}
}

void st_construct_unoptimised(input_data_t *p_input_data, alphabets_t *p_alphabets, st_node_t **pp_st_root){

	char *p_char_in_seq = p_input_data->pp_char_in_seq[0];
	int num_char_in_seq = p_input_data->p_num_char_in_seq[0];

	// Append $ at the end of sequence
	p_char_in_seq[num_char_in_seq] = '$';
	num_char_in_seq += 1;
	p_char_in_seq[num_char_in_seq] = 0;

	// check whether the last character is $
	if(p_char_in_seq[num_char_in_seq-1] != '$'){
		printf("Error is appending input string with $\n");
		exit(-1);
	}

	pp_st_root[0] = create_node(0, -1, (st_node_t *)0, (st_node_t *)0, 0, 0, 0, (st_node_t *)0);

	suffix_link_info_t s_sl_info; //dummy ignore for this function
	fill_suffix_link_info(&s_sl_info, pp_st_root[0], (st_node_t *)0);
	s_sl_info.p_v_node = (st_node_t *)0;
	s_sl_info.p_v_dash_node = (st_node_t *)0;
	int max_string_depth = 0;

	// For each suffix insert
	for(int i = 0; i < num_char_in_seq; i++){
		// Implement pushing
		find_path(pp_st_root[0], p_char_in_seq, i, num_char_in_seq, &s_sl_info, i);
		dfs(pp_st_root[0], p_input_data, &max_string_depth);
		__PRINTF__("End of loop %d\n", i);
	}

}

// Uses the U, U', V', beta node info to determine the V node and alpha for find_path function
void node_hop(suffix_link_info_t *p_sl_info, char *p_char_in_seq){
	st_node_t *p_u_node = p_sl_info->p_u_node;
	st_node_t *p_u_dash_node = p_sl_info->p_u_dash_node;

	// Checking if u_node is root
	if(p_u_node == p_u_node->p_suffix){
		if(p_u_node->node_idx != 0){
			printf("This should always happen for root node and not others \n");
			exit(-1);
		}
		p_sl_info->p_v_node = p_u_node;
		p_sl_info->p_v_dash_node = (st_node_t *)0;
		__PRINTF__("Node Hop Root %d\n", p_sl_info->p_v_node->node_idx);
		// p_sl_info->alpha = 0;
		return;
	}

	// Check if u_node has suffix link
	if(p_u_node->p_suffix){
		st_node_t *p_suffix = p_u_node->p_suffix;
		p_sl_info->p_v_node = p_suffix;
		p_sl_info->p_v_dash_node = (st_node_t *)0;
		// p_sl_info->alpha = p_suffix->edge_label[1] - p_suffix->edge_label[0];
		__PRINTF__("Node Hop U Node SL %d\n", p_sl_info->p_v_node->node_idx);
		return;
	}
	// Check if u_dash node has suffix link
	// If so have to do node hop
	if(p_u_dash_node->p_suffix){
		// If the u_dash is root node then offset by 1
		int u_start_idx = p_u_node->edge_label[0];
		if(p_u_dash_node->p_suffix == p_u_dash_node){
			if(p_u_dash_node->node_idx != 0){
				printf("This should always happen for root node and not others \n");
				exit(-1);
			}
			u_start_idx += 1;
		}
		st_node_t *p_suffix = p_u_dash_node->p_suffix;
		int beta = p_u_node->edge_label[1] - u_start_idx;
		// p_sl_info->alpha = 
		// parese through all the children to find the matching node
		st_node_t *p_prev_parent_node = (st_node_t *)0;
		st_node_t *p_parent_node = p_suffix;
		st_node_t *p_child_node = p_suffix->p_child;
		st_node_t *p_prev_sibling = (st_node_t *)0;
		// printf("Did it crash here???\n"); 
		while(beta){
			while(p_child_node){
				if(p_char_in_seq[p_child_node->edge_label[0]] == p_char_in_seq[u_start_idx]){
					int edge_length = p_child_node->edge_label[1] - p_child_node->edge_label[0];
					__PRINTF__("The beta value %d, %d, %d\n", beta, edge_length, p_child_node->node_idx);

					if(beta >= edge_length){
						for(int k=0; k<edge_length; k++){
							if(p_char_in_seq[p_child_node->edge_label[0]+k] != p_char_in_seq[u_start_idx+k]){
								printf("Error in matching %d, %c, %c\n", k, p_char_in_seq[p_child_node->edge_label[0]+k],
									p_char_in_seq[u_start_idx+k]);
								exit(-1);
							}
						}
						beta -= edge_length;
						u_start_idx += edge_length;
						p_prev_parent_node = p_parent_node;
						p_parent_node = p_child_node;
					}
					else{

						st_node_t *p_new_node = create_node(-1, -1,/*Leaf idx */ 
							(st_node_t *)0, /* No Suffix */
							p_child_node /*child */, 
							p_child_node->edge_label[0], 
							p_child_node->edge_label[0] + beta, 
							p_parent_node->string_depth, 
							p_child_node->p_sibling /*New node takes the child nodes siblings*/);

						if(p_prev_sibling){
							p_prev_sibling->p_sibling = p_new_node;
						}
						else{
							p_parent_node->p_child = p_new_node;
						}

						// Update the child node accordingly
						p_child_node->edge_label[0] += beta;
						p_child_node->p_sibling = (st_node_t *)0;

						// // create new node and make beta zero
						// // int parent_string_depth = p_child_node->string_depth - (p_child_node->edge_label[1] - p_child_node->edge_label[0]);
						// st_node_t *p_new_node = create_node(-1, p_child_node->leaf_idx, 
						// 	p_child_node->p_suffix, /* The lower node would store the suffix link information */
						// 	p_child_node->p_child /* child */, 
						// 	p_child_node->edge_label[0] + beta, 
						// 	p_child_node->edge_label[1], 
						// 	0 , // Pass dummy value which would be updated later , 
						// 	(st_node_t *)0 /* New node does not have any siblings */);
						// p_new_node->string_depth = p_child_node->string_depth;

						// // Update the child node accordingly
						// p_child_node->leaf_idx = -1;
						// p_child_node->p_child = p_new_node;
						// p_child_node->string_depth -= \
						// 	(p_child_node->edge_label[1]-p_child_node->edge_label[0]-beta);
						// p_child_node->edge_label[1] = p_child_node->edge_label[0] + beta;
						// p_child_node->p_suffix = (st_node_t *)0;

						// This is a recurrsive case where if the same character keeps happening again and again
						// we would have to split the node to make it its own suffix link.
						// When that happens the u node would shift to the lower node or the new node
						// While the upper node would be the suffix lik for it	
						// if(p_sl_info->p_u_node->node_idx == p_child_node->node_idx){
						// 	p_sl_info->p_u_node = p_new_node;
						// }

						// p_prev_parent_node = p_parent_node;
						// p_parent_node = p_child_node;// Update the child node to the new parent node
						p_prev_parent_node = p_parent_node;
						p_parent_node = p_new_node;
						beta = 0;
					}
					break;
				}

				if(p_child_node->p_sibling){
					p_prev_sibling = p_child_node;
					// iterate through all the siblings
					p_child_node = p_child_node->p_sibling;
				}
				// If the sibling does not exist create a node with the necessary sibling and insert
				else{
					printf("This condition should not happen %d\n", p_parent_node->node_idx);
					exit(-1);
					// This node is a intermediate node for new node to be added later
					// st_node_t *p_new_node = create_node(-1, -1, 
					// 		(st_node_t *)0, // No suffix info present
					// 		(st_node_t *)0,  // No child node present
					// 		u_start_idx, 
					// 		u_start_idx + beta, 
					// 		p_parent_node->string_depth, 
					// 		(st_node_t *)0 /* New node does not have any siblings */);
					// p_child_node->p_sibling = p_new_node;
					// p_child_node = p_child_node->p_sibling;
					// beta = 0;


					// p_prev_parent_node = p_parent_node;
					// p_parent_node = p_child_node;
					break;
				}
			}
			if(beta){
				p_prev_sibling = (st_node_t *)0;
				p_child_node = p_child_node->p_child;
			}
			else{
				break;
			}
		}
		//
		p_sl_info->p_v_dash_node = p_prev_parent_node;
		p_sl_info->p_v_node = p_parent_node;
		p_sl_info->p_u_node->p_suffix = p_parent_node;

		// if(p_sl_info->p_u_node->node_idx == p_sl_info->p_u_node->p_suffix->node_idx){
		// 	printf("How the hell can this happen???\n");
		// 	exit(-1);
		// }

		__PRINTF__("Node Hop U Node WO SL %d, %d\n", p_sl_info->p_v_node->node_idx, beta);
		// if(beta == 0){
		// 	p_u_node->p_suffix = p_parent_node;
		// }
		// else{
		// 	printf("Can this case happen ? %d\n", beta);
		// 	exit(-1);
		// }
		return;
	}

	printf("It should never come to this condition %d, %d\n", p_u_node->node_idx, p_u_dash_node->node_idx);
	exit(-1);
}

void st_construct(input_data_t *p_input_data, alphabets_t *p_alphabets, st_node_t **pp_st_root){

	char *p_char_in_seq = p_input_data->pp_char_in_seq[0];
	int num_char_in_seq = p_input_data->p_num_char_in_seq[0];

	// Append $ at the end of sequence
	p_char_in_seq[num_char_in_seq] = '$';
	num_char_in_seq += 1;
	p_char_in_seq[num_char_in_seq] = 0;
	p_input_data->p_num_char_in_seq[0] = num_char_in_seq;

	// check whether the last character is $
	if(p_char_in_seq[num_char_in_seq-1] != '$'){
		printf("Error is appending input string with $\n");
		exit(-1);
	}

	pp_st_root[0] = create_node(0, -1, (st_node_t *)0, (st_node_t *)0, 0, 0, 0, (st_node_t *)0);
	// For a root node the suffix link is same as the root node
	pp_st_root[0]->p_suffix = pp_st_root[0];
	
	// Initialise the suffix link info with root node
	suffix_link_info_t s_sl_info;
	fill_suffix_link_info(&s_sl_info, pp_st_root[0], (st_node_t *)0);
	s_sl_info.p_v_node = (st_node_t *)0;
	s_sl_info.p_v_dash_node = (st_node_t *)0;
	int max_string_depth = 0;

	clock_t begin = clock();

	// For each suffix do a node_hop to get to node v and then do a find path to insert a node
	for(int i = 0; i < num_char_in_seq; i++){
		for(int j=i; j<i+10; j++){
			__PRINTF__("%c, ", p_char_in_seq[j]);
		}
		__PRINTF__("\n");
		__PRINTF__("Current character index %d \n", i);
		// node_hop used to identify the node to start searching from
		node_hop(&s_sl_info, p_char_in_seq);
		// dfs(pp_st_root[0], p_char_in_seq);
		// printf("SL Info %x, %d\n", s_sl_info.p_v_node->node_idx, s_sl_info.p_v_node->string_depth+i);
		// Implement pushing
		find_path(s_sl_info.p_v_node, p_char_in_seq, s_sl_info.p_v_node->string_depth+i, 
			num_char_in_seq, &s_sl_info, i);
		// dfs(pp_st_root[0], p_char_in_seq);
		__PRINTF__("End of loop %d\n\n", i);
	}
	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

	printf("Total time for construction %f seconds\n", time_spent);	
	dfs(pp_st_root[0], p_input_data, &max_string_depth);
	printf("Maximum matching string %d\n", max_string_depth);
	printf("Total Number of Nodes %d, %ld, %ld, %ld\n", GLOBAL_NODE_IDX, sizeof(st_node_t), sizeof(int), sizeof(st_node_t*));
}