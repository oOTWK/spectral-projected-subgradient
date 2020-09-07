/*** 
Implementation of subgradient methods for set-covering problems (SCP)

1) spectral projected subgradient, based on
"Crema, A., Loreto, M., & Raydan, M. (2007) Spectral projected subgradient with a
momentum term for the Lagrangian dual approach. Computers and Operations Research,
34(10), 3174–3186."

2) basic subgradient, based on
"Beasley, J.E. (1990) A Lagrangian heuristic for set-covering problems. 
Naval Research Logistics, 37(1), 151–164."

***/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "subgradient.h"


#define FILE_FORMAT_ERR    fprintf(stderr, "Error: wrong SCP file format\n")
    
#define GETLINE(buf, buf_size, fp)  if (getline(&buf, &buf_size, fp) == -1) \
                                        { FILE_FORMAT_ERR; return -1; }

#define STR_TOKEN(token, s, delim)  if ((token = strtok(s, delim)) == NULL) \
                                        { FILE_FORMAT_ERR; return -1; }

#define MALLOC(var, type, size)     if ((var = (type) malloc(size)) == NULL) \
                                        { perror("Error malloc"); return -1; }

#define ZERO_TOL    pow(10, -12)


static int num_col, num_row, num_nonzero;
static int *costs;           // cost vector
static int *col_wise_a;      // column-wise constraint matrix
static int *col_wise_idx;    // start index of each column in colwise_a 
static int *row_wise_a;      // row-wise constraint matrix
static int *row_wise_idx;    // start index of each row in rowwise_a 
static int *col_sizes;
static int *row_sizes;
static double *best_dual_copy;


/* Initializes dual vector and computes its reduced cost and obj value.
Returns the initial obj value. */
static double init_dual_vector(double *dual, double *reduced_costs);

/* Computes subgradient vector (sps)
Returns -1 if current solution is optimal (i.e., subgradient vector becomes zero vector).
Returns 0 otherwise. */
static int compute_subg_vector_sps(int *subg, double *reduced_costs);

/* Computes subgradient vector (basic) 
Returns square norm of subgradient vector.
Returns -1 if current solution is optimal (i.e., subgradient vector becomes zero vector). */
static long long compute_subg_vector_basic(int *subg, double *reduced_costs, double *dual);



/* Reads SCP instance file and creates cost vector and constraint matrix.
Returns 0 on success, otherwise returns -1. */
int load_scp_instance(char *filename)
{
    int i, j, k;
    FILE *fp;
    char *buf = NULL;
    size_t buf_size = 0;
    char *token;

    if ((fp = fopen(filename, "r")) == NULL) { 
        perror("Error opening file"); return -1; 
    }

    // read first line: the number of row and the number of col
    GETLINE(buf, buf_size, fp)
    STR_TOKEN(token, buf, " ")
    num_row = atoi(token);
    STR_TOKEN(token, NULL, " ")
    num_col = atoi(token);

    if ((costs = (int *) malloc(num_col * sizeof(int))) == NULL) { 
        perror("Error malloc"); return -1;
    }

    // read cost vector
    for (i = 0, token = NULL; i < num_col; token = strtok(NULL, " ")) {
        if (token == NULL) {
            GETLINE(buf, buf_size, fp)
            STR_TOKEN(token, buf, " ")
        }

        if (*token != '\n') {
            costs[i++] = atoi(token);
        }
    }

    // temporary constraint matrix
    int **cols, **rows;
    MALLOC(cols, int **, num_col * sizeof(int *))
    MALLOC(rows, int **, num_row * sizeof(int *))
    for (i = 0; i < num_col; i++) {
        MALLOC(cols[i], int *, num_row * sizeof(int))
    }
    for (i = 0; i < num_row; i++) {
        MALLOC(rows[i], int *, num_col * sizeof(int))
    }

    MALLOC(col_sizes, int *, num_col * sizeof(int))
    memset(col_sizes, 0, num_col * sizeof(int));
    MALLOC(row_sizes, int *, num_row * sizeof(int))
    memset(row_sizes, 0, num_row * sizeof(int));

    // read rows (constraints)
    int col_idx;
    num_nonzero = 0;

    for (i = 0; i < num_row; i++) {
        GETLINE(buf, buf_size, fp)
        STR_TOKEN(token, buf, " ")
        row_sizes[i] = atoi(token);
        num_nonzero += row_sizes[i];

        for (j = 0, token = NULL; j < row_sizes[i]; token = strtok(NULL, " ")) {
            if (token == NULL) {
                GETLINE(buf, buf_size, fp)
                STR_TOKEN(token, buf, " ")
            }

            if (*token != '\n') {
                if ((col_idx = atoi(token) - 1) < 0) {
                    FILE_FORMAT_ERR; return -1;
                }
                cols[col_idx][col_sizes[col_idx]] = i;
                col_sizes[col_idx]++;
                rows[i][j++] = col_idx;
            }
        }
    }

    // create col-wise constraint matrix
    MALLOC(col_wise_a, int *, num_nonzero * sizeof(int))
    MALLOC(col_wise_idx, int *, (num_col+1) * sizeof(int))
    k = 0;
    for (i = 0; i < num_col; i++) {
        col_wise_idx[i] = k; // start index of i-th column
        for (j = 0; j < col_sizes[i]; j++) {
            col_wise_a[k++] = cols[i][j];
        }
    }
    col_wise_idx[num_col] = k;

    // create row-wise constraint matrix
    MALLOC(row_wise_a, int *, num_nonzero * sizeof(int))
    MALLOC(row_wise_idx, int *, (num_row+1) * sizeof(int))
    k = 0;
    for (i = 0; i < num_row; i++) {
        row_wise_idx[i] = k; // start index of i-th row
        for (j = 0; j < row_sizes[i]; j++) {
            row_wise_a[k++] = rows[i][j];
        }
    }
    row_wise_idx[num_row] = k;

    MALLOC(best_dual_copy, double *, num_row * sizeof(double))

    // free temporary constraint matrix
    for (i = 0; i < num_col; i++) {
        free(cols[i]);
    }
    for (i = 0; i < num_row; i++) {
        free(rows[i]);
    }
    free(cols);
    free(rows);

    return 0;
}


/* Initializes dual vector and computes its reduced cost and obj value.
Returns the initial obj value. */
static double init_dual_vector(double *dual, double *reduced_costs)
{
    int i, j, idx;
    double min_value, value, obj_value;

    obj_value = 0;

    // init dual_j = min (cost_i / size_i) for each col i in row j
    for (i = 0; i < num_row; i++) {
        min_value = costs[row_wise_a[row_wise_idx[i]]];
        for (j = row_wise_idx[i]; j < row_wise_idx[i+1]; j++) {
            idx = row_wise_a[j];
            value = (double) costs[idx] / col_sizes[idx];
            if (value < min_value) {
                min_value = value;
            }
        }
        dual[i] = min_value;
        obj_value += min_value;
    }

    // compute reduced cost
    for (i = 0; i < num_col; i++) {
        value = costs[i];
        for (j = col_wise_idx[i]; j < col_wise_idx[i+1]; j++) {
            idx = col_wise_a[j];
            value -= dual[idx];
        }
        reduced_costs[i] = value;
    }

    return obj_value;
}


/* Computes subgradient vector (sps)
Returns -1 if current solution is optimal (i.e., subgradient vector becomes zero vector).
Returns 0 otherwise. */
static int compute_subg_vector_sps(int *subg, double *reduced_costs)
{
    int i, j;

    for (i = 0; i < num_row; i++) {
        subg[i] = 1;
    }
    for (i = 0; i < num_col; i++) {
        if (reduced_costs[i] < pow(10,-14)) {
            for (j = col_wise_idx[i]; j < col_wise_idx[i+1]; j++) {
                subg[col_wise_a[j]]--;
            }
        }
    }

    for(i = 0; i < num_row; i++) {
        if (subg[i] != 0) {
            return 0;
        }
    }
    return -1;
}


/************** Spectral projected subgradient **************
Returns best (maximum) dual solution.
Returns -1 on system failure. */
double spectral_projected_subgradient(int max_itr)
{
    double curr_obj, best_obj, worst_obj, sub_obj, *past_objs;
    double *curr_dual, *old_dual, *best_dual, *dual1, *dual2;
    double *reduced_costs, *momentum, *dd;
    int worst_obj_idx, dd_size, *dd_idx;
    int *curr_subg, *old_subg, *subg1, *subg2;
    double alpha, alpha_deno, eta, eta_not, tau, accept, product, value;
    int itr, i, j, k;
    unsigned char is_opt;

    const int M = 10;
    const double mu = 0.7;
    const double gamma = 0.1;

    // allocate memory for local variables
    MALLOC(reduced_costs, double *, num_col * sizeof(double));
    MALLOC(dual1, double *, num_row * sizeof(double));
    MALLOC(dual2, double *, num_row * sizeof(double));
    MALLOC(subg1, int *, num_row * sizeof(int));
    MALLOC(subg2, int *, num_row * sizeof(int));

    MALLOC(past_objs, double *, M * sizeof(double));
    MALLOC(momentum, double *, num_row * sizeof(double));
    memset(momentum, 0, num_row * sizeof(double));
    MALLOC(dd, double *, num_row * sizeof(double));
    MALLOC(dd_idx, int *, num_row * sizeof(int)); // idx of nonzero values in vector dd
    
    // init data
    old_dual = curr_dual = dual1;
    best_dual = dual2;

    curr_obj = best_obj = worst_obj = past_objs[(worst_obj_idx=0)] 
    = init_dual_vector(curr_dual, reduced_costs);

    is_opt = compute_subg_vector_sps(subg1, reduced_costs);
    if (is_opt) goto cleanup;

    // compute eta_not
    eta_not = 0;
    for (i = 0; i < num_row; i++) {
        eta_not += subg1[i] * subg1[i];
    }
    eta_not = sqrt(eta_not);

    alpha = 0.1; // init alpha

    for (itr = 0; itr < max_itr; itr++) {
        // printf("%f\n", curr_obj);

        // swap subg
        if (itr % 2 == 0) {
            old_subg = subg1;
            curr_subg = subg2;
        } else {
            old_subg = subg2;
            curr_subg = subg1;
        }
        
        // update dual vector and objective value
        dd_size = 0;
        sub_obj = 0.0;
        product = 0.0;
        for (i = 0; i < num_row; i++) {
            momentum[i] = alpha * old_subg[i] + mu * momentum[i];
            value = old_dual[i] + momentum[i];
            if (value < 0) {
                value = 0;
            }
            value -= old_dual[i];
            if (value < - ZERO_TOL || value > ZERO_TOL) {
                dd[i]= value;
                product += value * momentum[i];

                curr_dual[i] = old_dual[i] + value;
                for (j = row_wise_idx[i]; j < row_wise_idx[i+1]; j++) {
                    reduced_costs[row_wise_a[j]] -= value;
                }

                dd_idx[dd_size] = i;
                dd_size++;
            } else {
                curr_dual[i] = old_dual[i];
            }
            sub_obj += curr_dual[i];
        }

        // compute current obj value
        curr_obj = sub_obj;
        for (i = 0; i < num_col; i++) {
            if (reduced_costs[i] < 0) {
                curr_obj += reduced_costs[i];
            }
        }

        // non-monotone line search along the direction dd
        product /= alpha;
        tau = 1.0;
        eta = eta_not / pow(itr, 1.1);
        accept = worst_obj + gamma * tau * product - eta;
        while (curr_obj < accept) {
            tau *= 0.5;
            // adjust dual vector
            for (k = 0; k < dd_size; k++) {
                i = dd_idx[k];
                value = tau * dd[i];
                if (value < - ZERO_TOL || value > ZERO_TOL) {
                    curr_dual[i] -= value;
                    sub_obj -= value;
                    for (j = row_wise_idx[i]; j < row_wise_idx[i+1]; j++) {
                        reduced_costs[row_wise_a[j]] += value;
                    }
                }
            }
            // compute adjusted obj value
            curr_obj = sub_obj;
            for (i = 0; i < num_col; i++) {
                if (reduced_costs[i] < 0) {
                    curr_obj += reduced_costs[i];
                }
            }
            accept -= gamma * tau * product;
        }

        // update best solution
        if (best_obj < curr_obj) {
            best_obj = curr_obj;

            // swap
            old_dual = curr_dual;
            curr_dual = best_dual;
            best_dual = old_dual;
        }  else {
            old_dual = curr_dual;
        }

        // compute subgradient vector
        is_opt = compute_subg_vector_sps(curr_subg, reduced_costs);
        if (is_opt) 
            break;

        // update alpha
        alpha = 0.0;
        alpha_deno = 0.0;
        for (k = 0; k < dd_size; k++) {
            i = dd_idx[k];
            value = dd[i];
            alpha += value * value;
            alpha_deno += value * (old_subg[i] - curr_subg[i]);
        }

        if (alpha_deno < ZERO_TOL) {
            alpha =  0.1;
        } else {
            alpha = tau * alpha / alpha_deno;
        }


        // update worst_lb
        i = (itr+1) % M;
        past_objs[i] = curr_obj;
        if (i == worst_obj_idx) {
            if (curr_obj <= worst_obj) {
                worst_obj = curr_obj;
            } else {
                worst_obj = curr_obj;
                for (j = M-1; j >= 0; j--) {
                    if (past_objs[j] < worst_obj) {
                        worst_obj = past_objs[j];
                        worst_obj_idx = j;
                    }
                }
            }
        } else if (curr_obj < worst_obj) {
            worst_obj = curr_obj;
            worst_obj_idx = i;
        }
    }


cleanup:
    if (is_opt) {
        best_dual = curr_dual;
        best_obj = curr_obj;
    }

    memcpy(best_dual_copy, best_dual, num_row * sizeof(double));

    free(reduced_costs);
    free(dual1);
    free(dual2);
    free(subg1);
    free(subg2);
    free(past_objs);
    free(momentum);
    free(dd);
    free(dd_idx);

    return best_obj;
}


/* Computes subgradient vector (basic) 
Returns square norm of subgradient vector.
Returns -1 if current solution is optimal.
(i.e., subgradient vector becomes zero vector) */
static long long compute_subg_vector_basic(int *subg, double *reduced_costs, double *dual)
{
    int i, j;
    long long norm;
    unsigned char is_opt;

    for (i = 0; i < num_row; i++) {
        subg[i] = 1;
    }
    for (i = 0; i < num_col; i++) {
        if (reduced_costs[i] < pow(10,-14)) {
            for (j = col_wise_idx[i]; j < col_wise_idx[i+1]; j++) {
                subg[col_wise_a[j]]--;
            }
        }
    }

    is_opt = 1;
    norm = 0;
    for(i = 0; i < num_row; i++) {
        if (subg[i] != 0) {
            is_opt = 0;

            if (subg[i] < 0 && dual[i] < pow(10,-14)) {
                subg[i] = 0;
            } else {
                norm += subg[i] * subg[i];
            }
        }
    }

    if (is_opt) {
        return -1;
    } else {
        if (norm == 0) norm = 1;
        return norm;
    }
}


/* Beasley's subgradient method 
Optimal upperbound (primal opt soln of original SCP) is given for test purpose.
Returns best (maximum) dual solution
Returns -1 on system failure */
double basic_subgradient(int max_itr, int upperbound)
{
    double curr_obj, best_obj;
    double *curr_dual, *old_dual, *best_dual, *dual1, *dual2;
    double *reduced_costs;
    int *subg;
    int itr, counter, i, j;
    long long norm;
    double lambda, step_size, value;

    const int counter_limit = 10;

    // allocate memory for local variables
    MALLOC(reduced_costs, double *, num_col * sizeof(double));
    MALLOC(dual1, double *, num_row * sizeof(double));
    MALLOC(dual2, double *, num_row * sizeof(double));
    MALLOC(subg, int *, num_row * sizeof(int));

    // init data
    old_dual = curr_dual = dual1;
    best_dual = dual2;

    curr_obj = best_obj = init_dual_vector(curr_dual, reduced_costs);

    itr = counter = 0;
    norm = 0;
    lambda = 2.0;
    for (itr = 0; itr < max_itr; itr++) {
        // printf("%f\n", curr_obj);

        // compute subgradient vector and step size
        norm = compute_subg_vector_basic(subg, reduced_costs, old_dual);
        if (norm < 0) 
            break;

        step_size = lambda * (1.05 * upperbound - curr_obj) / norm;

        // update dual vector and objective value
        curr_obj = 0.0;
        for (i = 0; i < num_row; i++) {
            value = step_size * subg[i];
            curr_dual[i] =  old_dual[i] + value;
            if (curr_dual[i] < 0) {
                value -= curr_dual[i];
                curr_dual[i] = 0.0;
            }

            if (value < - ZERO_TOL || value > ZERO_TOL) {
                for (j = row_wise_idx[i]; j < row_wise_idx[i+1]; j++) {
                    reduced_costs[row_wise_a[j]] -= value;
                }
            }
            curr_obj += curr_dual[i];
        }
        // compute current obj value
        for (i = 0; i < num_col; i++) {
            if (reduced_costs[i] < 0) {
                curr_obj += reduced_costs[i];
            }
        }

        // update best solution
        if (best_obj < curr_obj) {
            best_obj = curr_obj;
            counter = 0;

            // swap
            old_dual = curr_dual;
            curr_dual = best_dual;
            best_dual = old_dual;
        } else {
            counter++;
            old_dual = curr_dual;
        }

        if (counter > counter_limit) {
            lambda *= 0.5;
            counter = 0;
        }
    }


    if (norm < 0) {
        best_dual = curr_dual;
        best_obj = curr_obj;
    }
    
    memcpy(best_dual_copy, best_dual, num_row * sizeof(double));

    free(reduced_costs);
    free(dual1);
    free(dual2);
    free(subg);

    return best_obj;
}


// Copies best dual vector to the input dual.
void get_dual_vector(double *dual)
{
    memcpy(dual, best_dual_copy, num_row * sizeof(double));
}


// Computes reduced costs of best dual vector.
void get_reduced_costs(double *reduced_costs)
{
    int i, j;
    double value;
    for (i = 0; i < num_col; i++) {
        value = costs[i];
        for (j = col_wise_idx[i]; j < col_wise_idx[i+1]; j++) {
            value -= best_dual_copy[col_wise_a[j]];
        }
        reduced_costs[i] = value;
    }
}


int get_num_col() { return num_col; }
int get_num_row() { return num_row; }


void free_scp_instance()
{
    free(costs);
    free(col_wise_a);
    free(col_wise_idx);
    free(row_wise_a);
    free(row_wise_idx);
    free(col_sizes);
    free(row_sizes);
    free(best_dual_copy);
}
