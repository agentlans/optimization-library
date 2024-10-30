#include "opt.h"
#include <float.h>
#include <math.h>

#define SETERROR(err)                                                          \
  do {                                                                         \
    if (error)                                                                 \
      *error = (err);                                                          \
  } while (0)

// Function prototypes for forward declaration
static void swap(double *x, double *y);
static double compute_s(double a, double b, double c, double fa, double fb,
                        double fc);
static int should_bisect(double s, double a, double b, double c, double d,
                         int mflag);
static void update_interval(double *a, double *b, double *fa, double *fb,
                            double s, double fs);

double OPT_BrentRoot(double (*f)(double, void *), double a, double b,
                     void *params, double tol, int max_iter, OPT_Error *error) {

  // Compute initial function values
  double fa = f(a, params);
  double fb = f(b, params);

  // Check if the root is bracketed
  if (fa * fb >= 0) {
    SETERROR(OPT_ERROR_INVALID_BRACKET);
    return NAN; // Return NaN to indicate failure
  }

  // Ensure |f(a)| >= |f(b)|
  if (fabs(fa) < fabs(fb)) {
    swap(&a, &b);
    swap(&fa, &fb);
  }

  double c = a;   // Initialize c to a
  int mflag = 1;  // Flag for bisection method
  double d = 0.0; // Previous step size

  for (int iter = 0; iter < max_iter; iter++) {
    double fc = f(c, params); // Evaluate f at c

    // Inverse quadratic interpolation or secant method
    double s = compute_s(a, b, c, fa, fb, fc);

    // Check conditions for s and apply bisection if necessary
    if (should_bisect(s, a, b, c, d, mflag)) {
      s = (a + b) / 2; // Bisection method
      mflag = 1;       // Set mflag
    } else {
      mflag = 0; // Clear mflag
    }

    double fs = f(s, params); // Evaluate f at s

    // Update interval [a,b]
    update_interval(&a, &b, &fa, &fb, s, fs);

    // Ensure |f(a)| >= |f(b)| after each iteration
    if (fabs(fa) < fabs(fb)) {
      swap(&a, &b);
      swap(&fa, &fb);
    }

    // Check for convergence
    if ((fabs(b - a) < tol && fabs(fs) < tol) || fs == 0.0) {
      SETERROR(OPT_ERROR_SUCCESS);
      return s; // Return the root found within tolerance
    }
  }

  SETERROR(OPT_ERROR_MAX_ITERATIONS_REACHED); // Indicate that max iterations
                                              // were reached
  return NAN; // Return NaN to indicate failure to converge
}

// Helper function to swap two doubles
static void swap(double *x, double *y) {
  double temp = *x;
  *x = *y;
  *y = temp;
}

// Compute the next approximation using inverse quadratic interpolation
static double compute_inv_quad_s(double a, double b, double c, double fa,
                                 double fb, double fc) {
  double R = fb / fc;
  double S = fb / fa;
  double T = fa / fc;
  double P = S * (T * (R - T) * (c - b) - (1 - R) * (b - a));
  double Q = (T - 1) * (R - 1) * (S - 1);
  return b + P / Q;
}

static double compute_secant_s(double a, double b, double fa, double fb) {
  // Compute the next approximation using the secant method formula
  return b - fb * (b - a) / (fb - fa);
}

// Compute the next guess for the root using inverse quadratic interpolation or
// secant method
static double compute_s(double a, double b, double c, double fa, double fb,
                        double fc) {
  if (fa != fb && fa != fc && fb != fc) {
    return compute_inv_quad_s(a, b, c, fa, fb, fc);
  } else {
    return compute_secant_s(a, b, fa, fb);
  }
}

// Determine whether to use the bisection method based on various conditions
static int should_bisect(double s, double a, double b, double c, double d,
                         int mflag) {
  double lower_bound = (3 * a + b) / 4;
  double upper_bound = b;

  return (s < lower_bound || s > upper_bound ||
          (mflag && fabs(s - b) >= fabs(b - c) / 2) ||
          (!mflag && fabs(s - b) >= fabs(c - d) / 2));
}

// Update the interval [a,b] based on the value of s and fs
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
