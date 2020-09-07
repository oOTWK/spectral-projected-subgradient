/*** main.c for subgradient methods test on set-covering problems ***/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include "subgradient.h"

#define SPS 	1 // spectral projected subgradien
#define BASIC 	2 // basic subgradient

int main(int argc, char *argv[])
{	
	char *filename;
	clock_t begin_t, end_t;
	double dual_soln;
	int option;
	int upperbound;
	unsigned char subg_type = SPS;

	const int max_itr = 300;

	// parse option and get filename
	while ((option = getopt(argc, argv, "b:")) != -1) {
		if (option == 'b') {
			subg_type = BASIC;
			upperbound = atoi(optarg);
		} else {
			optind = argc;
			break;
		}
	}
	if (optind == argc) {
		fprintf(stderr, "usage: %s input_file [-b upperbound]\n", argv[0]);
		exit(1);
	}

	// read SCP file and init data structure
	filename = argv[optind];
	if (load_scp_instance(filename)) return 1;


	begin_t = clock();

	if (subg_type == SPS) {
		printf("Type: spectral projected subgradient\n");
		if ((dual_soln = spectral_projected_subgradient(max_itr)) < 0) return 1;
	} else {
		printf("Type: basic subgradient\n");
		if ((dual_soln = basic_subgradient(max_itr, upperbound)) < 0) return 1;
	}

	end_t = clock();

	printf("obj value: %f\n", dual_soln);
	printf("CPU time %.3f\n", (double) (end_t - begin_t) / CLOCKS_PER_SEC);


	/*** example: get dual vector and reduced costs ******
	int num_col, num_row;
	double *dual, *reduced_costs;
	num_col = get_num_col();
	num_row = get_num_row();
	dual = (double *) malloc(num_row * sizeof(double));
	reduced_costs = (double *) malloc(num_col * sizeof(double));
	get_dual_vector(dual);
	get_reduced_costs(reduced_costs);
	*********************************************************/


	free_scp_instance();

	return 0;
}