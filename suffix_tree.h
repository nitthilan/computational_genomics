#ifndef __SUFFIX_LINK_H__
#define __SUFFIX_LINK_H__

typedef struct _st_node_t_ {
	int node_idx;
	struct _st_node_t_ *p_suffix;
	struct _st_node_t_ *p_parent;
	int parent_label[2];
	int string_depth;
	struct _st_node_t_ *p_children;
}st_node_t;

void st_construct(input_data_t *p_input_data, alphabets_t *p_alphabets){

}
#endif //__SUFFIX_LINK_H__