#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "fasta_read.h"

// Read input sequence
void read_input_sequence(char *filename, input_data_t *p_input_data){

	FILE *fp;
	char * line = NULL;
    size_t len = 0;
    ssize_t read;
	fp = fopen(filename, "r");
	if(fp == NULL){
		printf("Input file not present. Could not open.\n");
		exit(-1);
	}

	int num_lines_in_seq[2];// num_char_in_seq[2], 
	// char *p_char_in_seq[2];
	int cur_seq_idx = -1, start_seq = 0;
	int MAX_CHAR_IN_LINE = 70;

	// Calculate the roung number of lines in the sequence
	// Since there are roughly 60 characters per line
	// this helps in allocating memory 
	while ((read = getline(&line, &len, fp)) != -1) {
		if(line[0] == '>'){
			cur_seq_idx += 1;
			num_lines_in_seq[cur_seq_idx] = 0;
			start_seq = 1;
			continue;
		}
		if(start_seq){
			num_lines_in_seq[cur_seq_idx] += 1;
		}
	}
	// Reset the variable and allocate memory for copying
	start_seq = 0; cur_seq_idx = -1;
	rewind(fp);
	for(int i = 0; i < 2; i++){
		p_input_data->p_char_in_seq[i] = (char *)malloc(MAX_CHAR_IN_LINE*num_lines_in_seq[i]*sizeof(char) + 2);
		if(p_input_data->p_char_in_seq[i] == NULL){
			printf("Cannot allocate memory for sequence (%d) copy. %d\n", i, num_lines_in_seq[i]);
			exit(-1);
		}
	}
	while ((read = getline(&line, &len, fp)) != -1) {
		// printf("%s, %lu, %ld\n", line, len, read);
		if(line[0] == '>'){
			cur_seq_idx += 1;
			p_input_data->num_char_in_seq[cur_seq_idx] = 0;
			start_seq = 1;
			continue;
		}
		if(start_seq){
			if(line[read-2] == '\r' && line[read-1] == '\n'){
				read -= 2; 
			}
			else if(line[read-1] == '\n'){
				read -= 1; 
			}
			while(read && line[read-1]==' '){
				read -= 1;
			}
			
			if(p_input_data->num_char_in_seq[cur_seq_idx]+read < \
				MAX_CHAR_IN_LINE*num_lines_in_seq[cur_seq_idx]){
				int num_char_in_seq = p_input_data->num_char_in_seq[cur_seq_idx];
				char *dst = &p_input_data->p_char_in_seq[cur_seq_idx][num_char_in_seq];
				char *src = line;
				memcpy(dst, src, read*sizeof(char));
				p_input_data->num_char_in_seq[cur_seq_idx] += read;
			}
			else{
				printf("Error in sequence read. %d, %d, %ld, %s\n",
					p_input_data->num_char_in_seq[cur_seq_idx], 
					num_lines_in_seq[cur_seq_idx], read,
					line);
				exit(-1);
			}
		}
        // printf("Retrieved line of length %zu:\n", read);
        // printf("%s", line);
    }

    // printf("List of sequences\n");
    // for(int i=0;i<2;i++){
    // 	p_char_in_seq[i][num_char_in_seq[i]] = '\0';
    // 	printf("%s\n",p_char_in_seq[i]);
    // }

	fclose(fp);
	if (line)
        free(line);

    return;
}

// read config file
void read_config(char *filename, param_config_t *p_param_config){
	FILE *fp;
	char * line = NULL;
    size_t len = 0;
    ssize_t read;
	fp = fopen(filename, "r");
	if(fp == NULL){
		printf("Input file not present. Could not open.\n");
		exit(-1);
	}

	char buf[100]; 
	int value;
    while (fscanf(fp, "%s %d", buf, &value)==2){
    	if(strcmp(buf, "match") == 0){
    		p_param_config->match = value;
    	}
    	if(strcmp(buf, "mismatch") == 0){
    		p_param_config->mismatch = value;
    	}
    	if(strcmp(buf, "g") == 0){
    		p_param_config->g = value;
    	}
    	if(strcmp(buf, "h") == 0){
    		p_param_config->h = value;
    	}
        // printf("%s, %d\n", buf, value);
    }
	fclose(fp);
	if (line)
        free(line);

    return;
}

void read_alphabets(char *filename, alphabets_t *p_alphabets){
	FILE *fp;
	char * line = NULL;
    size_t len = 0;
    ssize_t read;
	fp = fopen(filename, "r");
	if(fp == NULL){
		printf("Input file not present. Could not open.\n");
		exit(-1);
	}
	read = getline(&line, &len, fp);
	printf("%s, %lu, %ld\n", line, len, read);

	// Allocate memory to store the alphbets
	p_alphabets->p_char_in_seq = (char *)malloc(read*sizeof(char));
	if(p_alphabets->p_char_in_seq == NULL){
		printf("Cannot allocate memory for alphabets (%ld) copy\n", read);
		exit(-1);
	}
	for(int i = 0, j = 0; i < read; i++){
		if(line[i] != ' ' && line[i] != '\r' && line[i] != '\n'){
			p_alphabets->p_char_in_seq[j] = line[i];
			j+=1;
		}
		p_alphabets->num_alphabets = j;
	}
	printf("Num characters read %d\n", p_alphabets->num_alphabets);

	fclose(fp);
	if (line)
        free(line);
}

// Convert all symbols to small
void convert_to_small(input_data_t *p_input_data){
	for(int i=0; i<2; i++){
		for(int j=0; j<p_input_data->num_char_in_seq[i]; j++){
			// Value of "A" = 65. So if it is a capital letter
			// convert it to a small value
			if(p_input_data->p_char_in_seq[i][j] < (65+26)){
				p_input_data->p_char_in_seq[i][j] += 32;
			}
		}
	}
}


// Function to find the maximum between three values
max_t find_max(int a, int b, int c, int is_global){
	max_t max;
	if(a > b){
		if(a > c) {max.val = a; max.idx = 0;}
		else {max.val = c; max.idx = 2;}
	}
	else{
		if(b>c) {max.val = b; max.idx = 1;}
		else {max.val = c; max.idx = 2;}
	}
	if(max.val <= 0 && is_global == 0){
		max.val = 0;
		max.idx = -1;
	}
	return max;
}
