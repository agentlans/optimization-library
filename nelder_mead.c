#include "opt.h"
#include <math.h>
#include <stdlib.h>

// Macro for setting error code
#define SETERROR(err)                                                          \
  if (error)                                                                   \
  *error = (err)

// Constants for Nelder-Mead algorithm
const double ALPHA = 1.0; // Reflection coefficient
const double BETA = 0.5;  // Contraction coefficient
const double GAMMA = 2.0; // Expansion coefficient
const double DELTA = 0.5; // Shrink coefficient

// Function prototypes
static void calculate_centroid(double **simplex, double *centroid, int n);
static void order_simplex(double **simplex, double *function_values, int n);
static OPT_Error allocate_memory(double ***simplex, double **function_values,
                                 double **centroid, double **x_r, double **x_e,
                                 double **x_c, int n);
static void free_memory(double **simplex, double *function_values,
                        double *centroid, double *x_r, double *x_e, double *x_c,
                        int n);
static void update_simplex(double **simplex, double *function_values,
                           double *new_point, double new_value, int n);

// Main Nelder-Mead optimization function
void OPT_NelderMead(double (*f)(double *, int, void *), double *xmin, int n,
                    void *params, double tol, int max_iter, double initial_step,
                    OPT_Error *error) {
  double **simplex, *function_values, *centroid, *x_r, *x_e, *x_c;
  int iter_count = 0;

  // Allocate memory for simplex and other arrays
  OPT_Error err_code = allocate_memory(&simplex, &function_values, &centroid,
                                       &x_r, &x_e, &x_c, n);
  if (err_code != OPT_ERROR_SUCCESS) {
    SETERROR(err_code);
    return;
  }

  // Initialize simplex
  for (int i = 0; i <= n; i++) {
    for (int j = 0; j < n; j++) {
      simplex[i][j] = xmin[j] + (i == j ? initial_step : 0);
    }
    function_values[i] = f(simplex[i], n, params);
  }

  // Main optimization loop
  while (iter_count < max_iter) {
    // Order simplex vertices
    order_simplex(simplex, function_values, n);
    // Calculate centroid of n best points
    calculate_centroid(simplex, centroid, n);

    // Reflection
    for (int j = 0; j < n; j++) {
      x_r[j] = centroid[j] + ALPHA * (centroid[j] - simplex[n][j]);
    }
    double f_r = f(x_r, n, params);

    if (function_values[0] <= f_r && f_r < function_values[n - 1]) {
      update_simplex(simplex, function_values, x_r, f_r, n);
    } else if (f_r < function_values[0]) {
      // Expansion
      for (int j = 0; j < n; j++) {
        x_e[j] = centroid[j] + GAMMA * (x_r[j] - centroid[j]);
      }
      double f_e = f(x_e, n, params);
      update_simplex(simplex, function_values, f_e < f_r ? x_e : x_r,
                     f_e < f_r ? f_e : f_r, n);
    } else {
      // Contraction
      for (int j = 0; j < n; j++) {
        x_c[j] = centroid[j] + BETA * (simplex[n][j] - centroid[j]);
      }
      double f_c = f(x_c, n, params);
      if (f_c < function_values[n]) {
        update_simplex(simplex, function_values, x_c, f_c, n);
      } else {
        // Shrinkage
        for (int i = 1; i <= n; i++) {
          for (int j = 0; j < n; j++) {
            simplex[i][j] =
                simplex[0][j] + DELTA * (simplex[i][j] - simplex[0][j]);
          }
          function_values[i] = f(simplex[i], n, params);
        }
      }
    }

    // Check for convergence
    double variance_sum = 0.0, mean_value = 0.0;
    for (int i = 0; i <= n; i++) {
      mean_value += function_values[i];
    }
    mean_value /= (n + 1);

    for (int i = 0; i <= n; i++) {
      variance_sum += pow(function_values[i] - mean_value, 2);
    }

    if ((variance_sum / (n + 1)) < tol) {
      break;
    }

    iter_count++;
  }

  // Copy best solution to xmin
  for (int i = 0; i < n; i++) {
    xmin[i] = simplex[0][i];
  }

  // Free allocated memory
  free_memory(simplex, function_values, centroid, x_r, x_e, x_c, n);
  SETERROR((iter_count == max_iter) ? OPT_ERROR_MAX_ITERATIONS_REACHED
                                    : OPT_ERROR_SUCCESS);
}

// Update simplex with new point and function value
static void update_simplex(double **simplex, double *function_values,
                           double *new_point, double new_value, int n) {
  for (int j = 0; j < n; j++) {
    simplex[n][j] = new_point[j];
  }
  function_values[n] = new_value;
}

// Calculate centroid of the n best points
static void calculate_centroid(double **simplex, double *centroid, int n) {
  for (int j = 0; j < n; j++) {
    centroid[j] = 0.0;
    for (int i = 0; i < n; i++) {
      centroid[j] += simplex[i][j];
    }
    centroid[j] /= n;
  }
}

// Order simplex vertices based on function values
static void order_simplex(double **simplex, double *function_values, int n) {
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n + 1; j++) {
      if (function_values[i] > function_values[j]) {
        // Swap function values
        double temp_value = function_values[i];
        function_values[i] = function_values[j];
        function_values[j] = temp_value;

        // Swap simplex vertices
        double *temp_vertex = simplex[i];
        simplex[i] = simplex[j];
        simplex[j] = temp_vertex;
      }
    }
  }
}

// Allocate memory for simplex and other arrays
static OPT_Error allocate_memory(double ***simplex, double **function_values,
                                 double **centroid, double **x_r, double **x_e,
                                 double **x_c, int n) {
  *simplex = malloc((n + 1) * sizeof(double *));
  if (*simplex == NULL)
    return OPT_ERROR_MEMORY_ALLOCATION;

  for (int i = 0; i <= n; i++) {
    (*simplex)[i] = malloc(n * sizeof(double));
    if ((*simplex)[i] == NULL) {
      for (int j = 0; j < i; j++)
        free((*simplex)[j]);
      free(*simplex);
      return OPT_ERROR_MEMORY_ALLOCATION;
    }
  }

  *function_values = malloc((n + 1) * sizeof(double));
  *centroid = malloc(n * sizeof(double));
  *x_r = malloc(n * sizeof(double));
  *x_e = malloc(n * sizeof(double));
  *x_c = malloc(n * sizeof(double));

  if (*function_values == NULL || *centroid == NULL || *x_r == NULL ||
      *x_e == NULL || *x_c == NULL) {
    free_memory(*simplex, *function_values, *centroid, *x_r, *x_e, *x_c, n);
    return OPT_ERROR_MEMORY_ALLOCATION;
  }

  return OPT_ERROR_SUCCESS;
}

// Free allocated memory
static void free_memory(double **simplex, double *function_values,
                        double *centroid, double *x_r, double *x_e, double *x_c,
                        int n) {
  if (simplex) {
    for (int i = 0; i <= n; i++) {
      free(simplex[i]);
    }
    free(simplex);
  }
  free(function_values);
  free(centroid);
  free(x_r);
  free(x_e);
  free(x_c);
}
