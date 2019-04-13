#ifndef __SUFFIX_TREE_H__
#define __SUFFIX_TREE_H__

typedef struct _st_node_t_ {
	int node_idx; // Unique id for each node
	int leaf_idx; // Unique id for each suffix. All leaf less node had index -1
	struct _st_node_t_ *p_suffix;
	struct _st_node_t_ *p_child;
	int edge_label[2]; // start and end offset of the edge
	int string_depth; // Offset to the string from the root
	struct _st_node_t_ *p_sibling;
	int leaf_idx_range[2]; // Start and end idx of node
}st_node_t;

typedef struct _suffix_link_info_t_ {
	// Input for node_hop, Output for find path
	st_node_t *p_u_node;
	st_node_t *p_u_dash_node;
	// int beta; This info is present in u_node
	// Output for node_hop, Input for find path
	st_node_t *p_v_node;
	st_node_t *p_v_dash_node;
	// int alpha; This infor is present in v_node string_depth
}suffix_link_info_t;

void st_construct_unoptimised(input_data_t *p_input_data, alphabets_t *p_alphabets, st_node_t **pp_st_root);
void st_construct(input_data_t *p_input_data, alphabets_t *p_alphabets, st_node_t **pp_st_root);

#endif //__SUFFIX_TREE_H__