#pragma once
#include <cuda_runtime.h>
#include <math.h> /* log */

__device__ __forceinline__ double dot_product(const double *v1,
                                              const double *v2, int size) {
  double dot = 0.0;
  for (int i = 0; i < size; ++i) {
    dot += v1[i] * v2[i];
  }
  return dot;
}

template <unsigned int n>
__device__ __forceinline__ void householder(double *A, double *U) {
  for (int k = 0; k < n - 2; ++k) {
    double x[n]; // Adjust the size of x for necessary elements only

    // Copy the column below the diagonal to x
    for (int i = 0; i < n - k - 1; ++i) {
      x[i] = A[(k + 1 + i) * n + k];
    }

    // Compute norm of x, adjust x[0] for stability, and normalize x
    double x_norm = sqrt(dot_product(x, x, n - k - 1));
    x[0] = (x[0] > 0) ? (x[0] + x_norm) : (x[0] - x_norm);
    x_norm = sqrt(dot_product(x, x, n - k - 1));
    x_norm = (x_norm == 0) ? 1 : x_norm;

    for (int i = 0; i < n - k - 1; ++i) {
      x[i] /= x_norm;
    }

    // Apply Householder transformation H * A
    for (int i = 0; i < n; ++i) {
      double dot_product = 0.0;
      for (int j = k + 1; j < n; ++j) {
        dot_product += x[j - k - 1] * A[j * n + i]; // v^T A
      }
      for (int j = k + 1; j < n; ++j) {
        A[j * n + i] -= 2 * dot_product * x[j - k - 1]; // A - 2vv^T A
      }
    }

    // Apply Householder transformation A * H.T (since H is symmetric, H = H.T)
    for (int i = 0; i < n; ++i) {
      double dot_product = 0.0;
      for (int j = k + 1; j < n; ++j) {
        dot_product += A[i * n + j] * x[j - k - 1]; // A v
      }
      for (int j = k + 1; j < n; ++j) {
        A[i * n + j] -= 2 * dot_product * x[j - k - 1]; // A - 2 A v v^T
      }
    }

    // Now apply the transformation to U, which is just U * H
    for (int i = 0; i < n; ++i) {
      double dot_product = 0.0;
      for (int j = k + 1; j < n; ++j) {
        dot_product += U[i * n + j] * x[j - k - 1]; // U v
      }
      for (int j = k + 1; j < n; ++j) {
        U[i * n + j] -= 2 * dot_product * x[j - k - 1]; // U - 2 U v v^T
      }
    }
  }
}

template <unsigned int n>
__device__ __forceinline__ void householderTransformation(double *a,
                                                          double *u) {
  // General template, can be left empty or static_assert to cause a
  // compile-time error if this version gets instantiated with an unsupported
  // value of n
  householder<n>(a, u);
}
