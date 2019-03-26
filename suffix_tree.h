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
}st_node_t;

void st_construct(input_data_t *p_input_data, alphabets_t *p_alphabets, st_node_t **pp_st_root);

#endif //__SUFFIX_TREE_H__