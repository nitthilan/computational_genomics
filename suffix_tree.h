#ifndef __SUFFIX_LINK_H__
#define __SUFFIX_LINK_H__

typedef struct _st_node_t_ {
	int node_idx; // Index
	struct _st_node_t_ *p_suffix;
	struct _st_node_t_ *p_child;
	int edge_label[2];
	int string_depth;
	struct _st_node_t_ *p_sibling;
}st_node_t;

void st_construct(input_data_t *p_input_data, alphabets_t *p_alphabets, st_node_t *p_st_root);
#endif //__SUFFIX_LINK_H__