#ifndef __GENE_COMPARE_H__
#define __GENE_COMPARE_H__


// ./out ../Human-Mouse-BRCA2-cds.fasta 0 ../parameters.config
// Global Optimal Score 496
// Match 7686, Mismatch 2450, Gaps 1452, Opening gap 209 
// Match 7638, Mismatch 2496, Gaps 1456, Opening gap 137 
// Matching -8, -8, 0, 0



// Local score
// Local Optimal Score 1607
// Match 7377, Mismatch 2340, Gaps 512, Opening gap 129 
// Match 7351, Mismatch 2366, Gaps 512, Opening gap 100 

// ./out ../Opsin1_colorblindness_gene.fasta 0 ../parameters.config 
// Global Optimal Score -1471
// Match 2159, Mismatch 664, Gaps 1499, Opening gap 225 
// Match 2123, Mismatch 689, Gaps 1521, Opening gap 139 

// Precentage: Identities 0.499537, Gaps 0.346830
// Local
// Local Optimal Score 221
// Match 921, Mismatch 173, Gaps 210, Opening gap 33 
// Match 913, Mismatch 181, Gaps 210, Opening gap 24 

// Precentage: Identities 0.706288, Gaps 0.161043

// ./out ../sample.fasta 1 ../sample_parameters.config 
struct DP_cell {
    int *p_score;
    // ... // add any other field(s) that you may need for the implementation
};


#define S_IDX 0
#define I_IDX 2
#define D_IDX 1
typedef struct _cell_value_t {
	int cv[3]; // Value calculated for each direction
	char idx[3]; // The max idx of the prev cell
} cell_value_t;

typedef struct _dp_state{
	cell_value_t *p_dp_table; // Allocate 3xmxn to store the score for each cell
}dp_state_t;

// Storing output info for dumping
typedef struct _output_stat_t{
	int num_match;
	int num_mismatch;
	int num_gap;
	int num_open_gap;
	int opt_score;
	// int alignment_offset[2];
	char *p_output_trace; // Allocate 3x(m+n) to store the matching sequence
}output_stat_t;


void find_alignment(input_data_t *p_input_data, 
	dp_state_t *p_dp_state,
	param_config_t *p_param_config,
	int is_global,
	output_stat_t *p_output_stat);

void retrace_global(input_data_t *p_input_data,
	dp_state_t *p_dp_state, int is_global,
	int local_idx_i, int local_idx_j,
	output_stat_t *p_output_stat,
	param_config_t *p_param_config);

void print_output_stat(output_stat_t *p_output_stat);

#endif //__GENE_COMPARE_H__
