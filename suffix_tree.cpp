#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "fasta_read.h"
#include "suffix_tree.h"


typedef struct _st_node_t_ {
	int node_idx; // Unique id for each node
	int leaf_idx; // Unique id for each suffix. All leaf less node had index -1
	struct _st_node_t_ *p_suffix;
	struct _st_node_t_ *p_child;
	int edge_label[2];
	int string_depth;
	struct _st_node_t_ *p_sibling;
}st_node_t;

void find_path(st_node_t *p_node, char *p_char_in_seq, int offset){
	// Compare the edge label of each children
	p_sib_node = p_node;
	while(p_sib_node->p_sibling){

	}
}

// Wrapper function for initialising tha node
void init_node(st_node_t *p_root, int node_idx, int leaf_idx,
	st_node_t *p_suffix, st_node_t *p_child, int lbl_start_idx,
	int lbl_end_idx, int string_depth, st_node_t *p_sibling){
	// Initiliase node
	p_root->node_idx = node_idx;
	p_root->leaf_idx = leaf_idx;
	p_root->p_suffix = p_suffix;
	p_root->p_child = p_child;
	p_root->edge_label[0] = lbl_start_idx;
	p_root->edge_label[1] = lbl_end_idx;
	p_root->string_depth = string_depth;
	p_root->p_sibling = p_sibling;
}

void st_construct(input_data_t *p_input_data, alphabets_t *p_alphabets, st_node_t **pp_st_root){

	char *p_char_in_seq = p_input_data->p_char_in_seq[0];
	int num_char_in_seq = p_input_data->num_char_in_seq[0];

	// Append $ at the end of sequence
	p_char_in_seq[num_char_in_seq] = '$';
	num_char_in_seq += 1;

	// check whether the last character is $
	if(p_char_in_seq[num_char_in_seq-1] != '$'){
		printf("Error is appending input string with $\n");
		exit(-1);
	}

	// Add the full string as the first suffix
	pp_st_root[0] = (st_node_t *)malloc(1*sizeof(st_node_t));
	if(pp_st_root[0] == 0){
		printf("Error in allocating input node\n");
		exit(-1);
	}
	init_node(pp_st_root[0], 0, -1, 0, 0, 0, num_char_in_seq, 0, 0);


	// For each suffix insert
	for(int i = 0; i < num_char_in_seq; i++){
		// Implement pushing
		find_path(p_st_root, p_char_in_seq, i);
	}

}