
#ifndef __FASTA_READ_H__
#define __FASTA_READ_H__


// Storing input data
// typedef struct _input_data_t{
// 	char *p_char_in_seq[2];
// 	int num_char_in_seq[2];
// } input_data_t;
typedef struct _input_data_t{
	char **pp_char_in_seq;
	int *p_num_char_in_seq;
	int *p_num_lines_in_seq;
} input_data_t;
// Parameters used for weights
typedef struct _param_config_t{
	int match;
	int mismatch;
	int h;
	int g;
} param_config_t;

typedef struct _max_t {
	int val;
	int idx;
} max_t;
// Num alphabets
typedef struct _alphabets_t{
	char *p_char_in_seq;
	int num_alphabets;
}alphabets_t;

void convert_to_small(input_data_t *p_input_data);
max_t find_max(int a, int b, int c, int is_global);

void read_config(char *filename, param_config_t *p_param_config);
void read_input_sequence(char *filename, input_data_t *p_input_data, int num_input_to_skip,
	int num_seq_to_read);
void read_alphabets(char *filename, alphabets_t *p_alphabets);
int get_num_fasta(char *filename);

#endif // __FASTA_READ_H__