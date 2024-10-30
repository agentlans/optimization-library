#include "opt.h"
#include <math.h>
#include <stdio.h>

// Test function for OPT_BrentRoot
double test_function_root(double x, void *params) {
  return x * x - 4; // Root at x = 2 and x = -2
}

// Test function for OPT_BrentMinimize
double test_function_min(double x, void *params) {
  return (x - 2) * (x - 2) + 1; // Minimum at x = 2
}

// Test function for OPT_NelderMead
double test_function_nelder_mead(double *x, int n, void *params) {
  return 100 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) +
         (1 - x[0]) * (1 - x[0]); // Rosenbrock function
}

int main() {
  OPT_Error error;

  // Test OPT_BrentRoot
  printf("Testing OPT_BrentRoot:\n");
  double root =
      OPT_BrentRoot(test_function_root, 0, 3, NULL, 1e-6, 100, &error);
  if (error == OPT_ERROR_SUCCESS) {
    printf("Root found: %f\n", root);
  } else {
    printf("Error in OPT_BrentRoot: %d\n", error);
  }

  // Test OPT_BrentMinimize
  printf("\nTesting OPT_BrentMinimize:\n");
  double xmin;
  double min_value = OPT_BrentMinimize(test_function_min, 0, 4, NULL, 1e-6, 100,
                                       &xmin, &error);
  if (error == OPT_ERROR_SUCCESS) {
    printf("Minimum found at x = %f with value %f\n", xmin, min_value);
  } else {
    printf("Error in OPT_BrentMinimize: %d\n", error);
  }

  // Test OPT_NelderMead
  printf("\nTesting OPT_NelderMead:\n");
  double x[2] = {-1.2, 1.0}; // Initial guess
  OPT_NelderMead(test_function_nelder_mead, x, 2, NULL, 1e-6, 1000, 1.0,
                 &error);
  if (error == OPT_ERROR_SUCCESS) {
    printf("Minimum found at (%f, %f)\n", x[0], x[1]);
    printf("Function value at minimum: %f\n",
           test_function_nelder_mead(x, 2, NULL));
  } else {
    printf("Error in OPT_NelderMead: %d\n", error);
  }

  return 0;
}
