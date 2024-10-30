#include <float.h>
#include <math.h>
#include "opt.h"

static const double CGOLD = 0.3819660;
static const double ZEPS = 1.0e-10;
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define SETERROR(err)                                                          \
  do {                                                                         \
    if (error)                                                                 \
      *error = (err);                                                          \
  } while (0)

double OPT_BrentMinimize(double (*f)(double, void *), double a, double b,
                         void *params, double tol, int max_iter, double *xmin,
                         OPT_Error *error) {
  int iter;
  double x, w, v, u, fx, fw, fv, fu;
  double q, r, p, etemp, d, e;
  double xm, tol1, tol2;
  etemp = 0;
  e = 0;

  SETERROR(OPT_ERROR_SUCCESS);

  // Ensure a < b
  if (a > b) {
    double tmp = a;
    a = b;
    b = tmp;
  }

  // Initialize x, w, and v
  x = a + CGOLD * (b - a); // Start with golden section
  w = v = x;
  fx = fw = fv = (*f)(x, params);

  // Check endpoints
  double fa = (*f)(a, params);
  double fb = (*f)(b, params);
  if (fa < fx) {
    x = a;
    fx = fa;
  }
  if (fb < fx) {
    x = b;
    fx = fb;
  }

  for (iter = 1; iter <= max_iter; iter++) {
    xm = 0.5 * (a + b);
    tol1 = tol * fabs(x) + ZEPS;
    tol2 = 2.0 * tol1;

    // Check if done
    if (fabs(x - xm) <= (tol2 - 0.5 * (b - a))) {
      *xmin = x;
      return fx;
    }

    if (fabs(e) > tol1) {
      // Parabolic interpolation
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = 2.0 * (q - r);
      if (q > 0.0)
        p = -p;
      q = fabs(q);
      etemp = e;
      e = d;

      if (fabs(p) < fabs(0.5 * q * etemp) && p > q * (a - x) &&
          p < q * (b - x)) {
        // Parabolic step
        d = p / q;
        u = x + d;
        if (u - a < tol2 || b - u < tol2)
          d = SIGN(tol1, xm - x);
      } else {
        // Golden section step
        e = (x >= xm) ? a - x : b - x;
        d = CGOLD * e;
      }
    } else {
      // Golden section step
      e = (x >= xm) ? a - x : b - x;
      d = CGOLD * e;
    }

    // Update u
    u = (fabs(d) >= tol1) ? x + d : x + SIGN(tol1, d);
    fu = (*f)(u, params);

    // Update a, b, v, w, and x
    if (fu <= fx) {
      if (u >= x)
        a = x;
      else
        b = x;
      v = w; // Track the previous best
      w = x;
      x = u;
      fv = fw; // Update the best function values
      fw = fx;
      fx = fu;
    } else {
      if (u < x)
        a = u;
      else
        b = u;

      // Ensure proper tracking of the previous best values
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

  SETERROR(OPT_ERROR_MAX_ITERATIONS_REACHED);
  *xmin = x;
  return fx;
}
