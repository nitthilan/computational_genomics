#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "fasta_read.h"
#include "suffix_tree.h"


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


int GLOBAL_NODE_IDX = 0;
// Wrapper function for initialising tha node
st_node_t *create_node(int node_idx, int leaf_idx,
	st_node_t *p_suffix, st_node_t *p_child, int lbl_start_idx,
	int lbl_end_idx, int string_depth, st_node_t *p_sibling){

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
	p_node->string_depth = string_depth;
	p_node->p_sibling = p_sibling;
	return p_node;
}

void find_path(st_node_t *p_root_node, char *p_char_in_seq, int leaf_offset, int size_of_seq){

	st_node_t *p_parent_node = p_root_node;
	int cur_char_offset = leaf_offset;
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
			st_node_t *p_new_node = create_node(-1, leaf_offset, (st_node_t *)0 /*TBD*/, 
				(st_node_t *)0 /*child */, 
				cur_char_offset, size_of_seq, p_parent_node->string_depth+1, 
				p_sib_node);
			if(p_prev_sib_node){
				// if p_prev_sib_node is not null then add the current node as the next sibling
				p_prev_sib_node->p_sibling = p_new_node;
			}
			else{
				// if p_prev_sib_node is null then add the current node as the leftmost node
				// Make the new node the first sibling
				p_parent_node->p_child = p_new_node;
			}
			break;
		}
		// is_match_found = 2 : Its a partial match. So break and introduce a new node
		if(is_match_found == 2){
			// Create a intemediate node and introduce the sequence as its child
			st_node_t *p_new_node = create_node(-1, p_sib_node->leaf_idx, (st_node_t *)0 /*TBD*/, 
				p_sib_node->p_child /*child */, 
				p_sib_node->edge_label[0] + partial_match_offset, 
				p_sib_node->edge_label[1], p_sib_node->string_depth+1, 
				(st_node_t *)0 /*Update this after comparing between the two new nodes */);
			// Update the sibling with the partial node information
			p_sib_node->edge_label[1] = p_sib_node->edge_label[0] + partial_match_offset;
			p_sib_node->leaf_idx = -1;

			st_node_t *p_new_leaf_node = create_node(-1, leaf_offset, (st_node_t *)0 /*TBD*/, 
				(st_node_t *)0 /*child */, 
				cur_char_offset + partial_match_offset, size_of_seq, p_sib_node->string_depth+1, 
				(st_node_t *)0 /*Update this after comparing between the two new nodes */);
			// Maintaining the ordering of siblings
			if(p_char_in_seq[p_sib_node->edge_label[0] + partial_match_offset] < 
			   p_char_in_seq[cur_char_offset + partial_match_offset]){
			   	p_sib_node->p_child = p_new_node;
				p_new_node->p_sibling = p_new_leaf_node;
				p_new_leaf_node->p_sibling = (st_node_t *)0;
			}
			else{
				p_sib_node->p_child = p_new_leaf_node;
				p_new_leaf_node->p_sibling = p_new_node;
				p_new_node->p_sibling = (st_node_t *)0;
			}
			break;
		}

		// is_match_found = 1 : Its a perfect match and so need to goto the next node
		if(is_match_found == 1){
			// Make the current sibling as the next parent and go down the tree
			p_parent_node = p_sib_node;
			cur_char_offset += (end_idx-start_idx);
		}
	}
}

void st_construct(input_data_t *p_input_data, alphabets_t *p_alphabets, st_node_t **pp_st_root){

	char *p_char_in_seq = p_input_data->p_char_in_seq[0];
	int num_char_in_seq = p_input_data->num_char_in_seq[0];

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


	// For each suffix insert
	for(int i = 0; i < num_char_in_seq; i++){
		// Implement pushing
		find_path(pp_st_root[0], p_char_in_seq, i, num_char_in_seq);
	}

}