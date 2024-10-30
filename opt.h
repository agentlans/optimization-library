#ifndef OPT_H
#define OPT_H

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  OPT_ERROR_SUCCESS = 0,         /**< Indicates successful execution. */
  OPT_ERROR_INVALID_BRACKET = 1, /**< Indicates that the root is not bracketed
                                    or the bracket is invalid. */
  OPT_ERROR_MAX_ITERATIONS_REACHED =
      2, /**< Indicates that the maximum number of iterations was reached. */
  OPT_ERROR_MEMORY_ALLOCATION =
      3,                       /**< Indicates a memory allocation failure. */
  OPT_ERROR_INVALID_INPUT = 4, /**< Indicates invalid input parameters. */
  OPT_ERROR_UNKNOWN = 5 /**< Indicates an unknown or unspecified error. */
} OPT_Error;

/**
 * Finds the root of a function using Brent's method.
 *
 * @param f The function to find the root of.
 * @param x1 The lower bound of the initial bracket.
 * @param x2 The upper bound of the initial bracket.
 * @param params Additional parameters to pass to the function.
 * @param tol The desired tolerance for the root.
 * @param max_iter The maximum number of iterations allowed.
 * @param error Pointer to store the error code.
 * @return The estimated root of the function.
 */
double OPT_BrentRoot(double (*f)(double, void *), double x1, double x2,
                     void *params, double tol, int max_iter, OPT_Error *error);

/**
 * Minimizes a function using Brent's method.
 *
 * @param f The function to minimize.
 * @param a The lower bound of the search interval.
 * @param b The upper bound of the search interval.
 * @param params Additional parameters to pass to the function.
 * @param tol The desired tolerance for the minimum.
 * @param max_iter The maximum number of iterations allowed.
 * @param xmin Pointer to store the x-value at the minimum.
 * @param error Pointer to store the error code.
 * @return The minimum value of the function.
 */
double OPT_BrentMinimize(double (*f)(double, void *), double a, double b,
                         void *params, double tol, int max_iter, double *xmin,
                         OPT_Error *error);

/**
 * Minimizes a multidimensional function using the Nelder-Mead method.
 *
 * @param f The function to minimize.
 * @param xmin Array to store the coordinates of the minimum.
 * @param n The number of dimensions.
 * @param params Additional parameters to pass to the function.
 * @param tol The desired tolerance for the minimum.
 * @param max_iter The maximum number of iterations allowed.
 * @param initial_step The initial step size for creating the simplex.
 * @param error Pointer to store the error code.
 */
void OPT_NelderMead(double (*f)(double *, int, void *), double *xmin, int n,
                    void *params, double tol, int max_iter, double initial_step,
                    OPT_Error *error);

#ifdef __cplusplus
}
#endif

#endif /* OPT_H */
