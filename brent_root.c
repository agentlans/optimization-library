#include "opt.h"
#include <float.h>
#include <math.h>

// Macro to set error if error pointer is not NULL
#define SETERROR(err)                                                          \
  if (error)                                                                   \
  *error = (err)

// Function prototypes
static void swap(double *x, double *y);
static double compute_s(double a, double b, double c, double fa, double fb,
                        double fc);
static int should_bisect(double s, double a, double b, double c, double d,
                         int mflag);
static void update_interval(double *a, double *b, double *fa, double *fb,
                            double s, double fs);

// Brent's method for root finding
double OPT_BrentRoot(double (*f)(double, void *), double a, double b,
                     void *params, double tol, int max_iter, OPT_Error *error) {
  // Compute function values at interval endpoints
  double fa = f(a, params), fb = f(b, params);

  // Check if root is bracketed
  if (fa * fb >= 0) {
    SETERROR(OPT_ERROR_INVALID_BRACKET);
    return NAN;
  }

  // Ensure |f(a)| >= |f(b)|
  if (fabs(fa) < fabs(fb)) {
    swap(&a, &b);
    swap(&fa, &fb);
  }

  double c = a, d = 0.0, fc, fs, s;
  int mflag = 1; // Flag for bisection method

  for (int iter = 0; iter < max_iter; iter++) {
    fc = f(c, params);

    // Compute next approximation
    s = compute_s(a, b, c, fa, fb, fc);

    // Check if bisection is needed
    if (should_bisect(s, a, b, c, d, mflag)) {
      s = (a + b) / 2; // Use bisection
      mflag = 1;
    } else {
      mflag = 0;
    }

    fs = f(s, params);
    d = c;
    c = b;

    // Update interval [a, b]
    update_interval(&a, &b, &fa, &fb, s, fs);

    // Ensure |f(a)| >= |f(b)|
    if (fabs(fa) < fabs(fb)) {
      swap(&a, &b);
      swap(&fa, &fb);
    }

    // Check for convergence
    if ((fabs(b - a) < tol && fabs(fs) < tol) || fs == 0.0) {
      SETERROR(OPT_ERROR_SUCCESS);
      return s;
    }
  }

  SETERROR(OPT_ERROR_MAX_ITERATIONS_REACHED);
  return NAN;
}

// Swap two double values
static void swap(double *x, double *y) {
  double temp = *x;
  *x = *y;
  *y = temp;
}

// Compute next approximation using inverse quadratic interpolation or secant
// method
static double compute_s(double a, double b, double c, double fa, double fb,
                        double fc) {
  if (fa != fb && fa != fc && fb != fc) {
    // Inverse quadratic interpolation
    double R = fb / fc, S = fb / fa, T = fa / fc;
    double P = S * (T * (R - T) * (c - b) - (1 - R) * (b - a));
    double Q = (T - 1) * (R - 1) * (S - 1);
    return b + P / Q;
  } else {
    // Secant method
    return b - fb * (b - a) / (fb - fa);
  }
}

// Determine if bisection should be used
static int should_bisect(double s, double a, double b, double c, double d,
                         int mflag) {
  double lower_bound = (3 * a + b) / 4;
  return (s < lower_bound || s > b ||
          (mflag && fabs(s - b) >= fabs(b - c) / 2) ||
          (!mflag && fabs(s - b) >= fabs(c - d) / 2));
}

// Update the interval [a, b] based on new approximation
static void update_interval(double *a, double *b, double *fa, double *fb,
                            double s, double fs) {
  if (*fa * fs < 0) {
    *b = s;
    *fb = fs;
  } else {
    *a = s;
    *fa = fs;
  }
}
