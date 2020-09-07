/*** header file, subgradient methods for set-covering problems

1) spectral projected subgradient, based on
"Crema, A., Loreto, M., & Raydan, M. (2007) Spectral projected subgradient with a
momentum term for the lagrangian dual approach. Computers and OperationsResearch,
34(10), 3174–3186."

2) basic subgradient, based on
"Beasley, J.E. (1990) A Lagrangian heuristic for set-covering problems. 
Naval Research Logistics, 37(1), 151–164."

***/

#ifndef Subgradient_h
#define Subgradient_h

/* Reads SCP instance file and creates cost vector and constraint matrix.
Returns 0 on success, otherwise returns -1. */
int load_scp_instance(char *filename);

/* Spectral projected subgradient 
Returns best (maximum) dual solution.
Returns -1 on system failure. */
double spectral_projected_subgradient(int max_itr);

/* Beasley's subgradient method 
Optimal upperbound (primal opt soln of original SCP) is given for test purpose.
Returns best (maximum) dual solution
Returns -1 on system failure */
double basic_subgradient(int max_itr, int upperbound);

// Copies best dual vector to the input dual.
void get_dual_vector(double *dual);

// Computes reduced costs of best dual vector.
void get_reduced_costs(double *reduced_costs);

int get_num_col();
int get_num_row();

void free_scp_instance();

#endif /* Subgradient_h */
