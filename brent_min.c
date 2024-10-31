#include "opt.h"
#include <float.h>
#include <math.h>

#define CGOLD 0.3819660 // Golden ratio constant
#define ZEPS 1.0e-10    // Small number to prevent division by zero
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define SETERROR(err)                                                          \
  if (error)                                                                   \
  *error = (err)

double OPT_BrentMinimize(double (*f)(double, void *), double a, double b,
                         void *params, double tol, int max_iter, double *xmin,
                         OPT_Error *error) {
  double x, w, v, fx, fw, fv, u, fu;
  double e = 0.0, d = 0.0;

  SETERROR(OPT_ERROR_SUCCESS);

  // Ensure a < b
  if (a > b) {
    double tmp = a;
    a = b;
    b = tmp;
  }

  // Initialize x, w, and v using golden section
  x = w = v = a + CGOLD * (b - a);
  fx = fw = fv = f(x, params);

  // Check if endpoints are better
  double fa = f(a, params), fb = f(b, params);
  if (fa < fx) {
    x = a;
    fx = fa;
  }
  if (fb < fx) {
    x = b;
    fx = fb;
  }

  for (int iter = 1; iter <= max_iter; iter++) {
    double xm = 0.5 * (a + b);
    double tol1 = tol * fabs(x) + ZEPS;
    double tol2 = 2.0 * tol1;

    // Check if we've reached the desired tolerance
    if (fabs(x - xm) <= (tol2 - 0.5 * (b - a))) {
      *xmin = x;
      return fx;
    }

    // Attempt parabolic fit
    if (fabs(e) > tol1) {
      // Compute parabolic interpolation
      double r = (x - w) * (fx - fv);
      double q = (x - v) * (fx - fw);
      double p = (x - v) * q - (x - w) * r;
      q = 2.0 * (q - r);
      if (q > 0.0)
        p = -p;
      q = fabs(q);
      double etemp = e;
      e = d;

      // Check if parabolic fit is acceptable
      if (fabs(p) < fabs(0.5 * q * etemp) && p > q * (a - x) &&
          p < q * (b - x)) {
        d = p / q;
        u = x + d;
        // Make sure u is not too close to a or b
        if (u - a < tol2 || b - u < tol2)
          d = SIGN(tol1, xm - x);
      } else {
        // Fall back to golden section
        e = (x >= xm) ? a - x : b - x;
        d = CGOLD * e;
      }
    } else {
      // Use golden section
      e = (x >= xm) ? a - x : b - x;
      d = CGOLD * e;
    }

    // Update u, ensuring it's at least tol1 away from x
    u = (fabs(d) >= tol1) ? x + d : x + SIGN(tol1, d);
    fu = f(u, params);

    // Update the interval and best point
    if (fu <= fx) {
      if (u >= x)
        a = x;
      else
        b = x;
      v = w;
      w = x;
      x = u;
      fv = fw;
      fw = fx;
      fx = fu;
    } else {
      if (u < x)
        a = u;
      else
        b = u;
      if (fu <= fw || w == x) {
        v = w;
        w = u;
        fv = fw;
        fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v = u;
        fv = fu;
      }
    }
  }

  // Maximum iterations reached
  SETERROR(OPT_ERROR_MAX_ITERATIONS_REACHED);
  *xmin = x;
  return fx;
}
