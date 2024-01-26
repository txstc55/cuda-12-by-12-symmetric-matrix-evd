#pragma once
#include "householder.cuh"
#include "qr_tri_diagonal.cuh"
#include <cuda_runtime.h>
// this is an attempt to solve for eigen vectors directly from the tri diagonal
// matrix it's working and not working yet the problem is the large error it got
// from the substitution
template <unsigned int n>
__device__ __forceinline__ void
computeEigenVectors(const double A[n * n], const double v[n], double *ev) {
  double AForward[n * n];
  double ABackward[n * n];
  double chosenA[n * n];
  for (int i = 0; i < n; i++) {
    double minimumError = 0.0;
    bool forward = true;
    double lambda = v[i];
    for (unsigned int j = 0; j < n * n; j++) {
      AForward[j] = A[j];
    }
    for (int j = 0; j < n; j++) {
      AForward[j * n + j] -= lambda;
    }

    // we will implement a row swapping gauss elimination
    {
      // do the gauss elimination
      double factor = AForward[n] / AForward[0];
      for (int j = 0; j < n; j++) {
        AForward[j + n] -= factor * AForward[j];
      }
    }

    for (int j = 1; j < n - 1; j++) {
      // now perform the gauss elimination
      double factor = AForward[(j + 1) * n + j] / AForward[j * n + j];
      for (int k = j + 1; k < n; k++) {
        AForward[(j + 1) * n + k] -= factor * AForward[j * n + k];
      }
      AForward[(j + 1) * n + j] = 0.0;
    }
    minimumError = fabs(AForward[n * n - 1]);

    // now try the backward version
    for (unsigned int j = 0; j < n * n; j++) {
      ABackward[j] = A[n * n - 1 - j];
    }
    for (int j = 0; j < n; j++) {
      ABackward[j * n + j] -= lambda;
    }
    // we will implement a row swapping gauss elimination
    {
      // do the gauss elimination
      double factor = ABackward[n] / ABackward[0];
      for (int j = 0; j < n; j++) {
        ABackward[j + n] -= factor * ABackward[j];
      }
    }

    for (int j = 1; j < n - 1; j++) {
      // now perform the gauss elimination
      double factor = ABackward[(j + 1) * n + j] / ABackward[j * n + j];
      for (int k = j + 1; k < n; k++) {
        ABackward[(j + 1) * n + k] -= factor * ABackward[j * n + k];
      }
      ABackward[(j + 1) * n + j] = 0.0;
    }

    unsigned int chosenPositions[n];
    if (fabs(ABackward[n * n - 1]) < minimumError) {
      minimumError = fabs(ABackward[n * n - 1]);
      for (unsigned int j = 0; j < n * n; j++) {
        chosenA[j] = ABackward[j];
      }
      for (unsigned int j = 0; j < n; j++) {
        chosenPositions[j] = n - 1 - j;
      }
      forward = false;
    } else {
      for (unsigned int j = 0; j < n * n; j++) {
        chosenA[j] = AForward[j];
      }
      for (unsigned int j = 0; j < n; j++) {
        chosenPositions[j] = j;
      }
    }

    // ok actually do the back substitution
    double solution[n];
    solution[n - 1] = 1.0;
    for (int j = n - 2; j >= 0; j--) {
      double sum = 0.0;
      for (unsigned int k = j + 1; k < n; k++) {
        sum += chosenA[j * n + k] * solution[k];
      }
      solution[j] = -sum / chosenA[j * n + j];
    }
    // normalize the solution
    double sum = 0.0;
    for (unsigned int j = 0; j < n; j++) {
      sum += solution[j] * solution[j];
    }
    sum = sqrt(sum);
    for (unsigned int j = 0; j < n; j++) {
      solution[j] /= sum;
    }

    // here we do 2 iterations of power iteration
    for (int iteration = 0; iteration < 2; iteration++) {
      double newSolution[n];
      // first copy the matrix
      if (!forward) {
        for (unsigned int j = 0; j < n * n; j++) {
          chosenA[j] = A[n * n - 1 - j];
        }
        for (int j = 0; j < n; j++) {
          chosenA[j * n + j] -= lambda;
        }
      } else {
        for (unsigned int j = 0; j < n * n; j++) {
          chosenA[j] = A[j];
        }
        for (int j = 0; j < n; j++) {
          chosenA[j * n + j] -= lambda;
        }
      }

      // now do the row swapping gauss elimination
      {
        // deal with the first row
        if (fabs(chosenA[0]) < fabs(chosenA[n])) {
          // swap the first row with the second row
          for (int j = 0; j < n; j++) {
            double temp = chosenA[j];
            chosenA[j] = chosenA[j + n];
            chosenA[j + n] = temp;
          }
          double temp = solution[0];
          solution[0] = solution[1];
          solution[1] = temp;
        }
        // do the gauss elimination
        double factor = chosenA[n] / chosenA[0];
        for (int j = 0; j < n; j++) {
          chosenA[j + n] -= factor * chosenA[j];
        }
        solution[1] -= factor * solution[0];
      }

      for (int j = 1; j < n - 1; j++) {
        if (fabs(chosenA[j * n + j]) < fabs(chosenA[(j + 1) * n + j])) {
          // swap the rows
          for (int k = 0; k < n; k++) {
            double temp = chosenA[j * n + k];
            chosenA[j * n + k] = chosenA[(j + 1) * n + k];
            chosenA[(j + 1) * n + k] = temp;
          }
          double temp = solution[j];
          solution[j] = solution[j + 1];
          solution[j + 1] = temp;
        }
        // now perform the gauss elimination
        double factor = chosenA[(j + 1) * n + j] / chosenA[j * n + j];
        for (int k = j + 1; k < n; k++) {
          chosenA[(j + 1) * n + k] -= factor * chosenA[j * n + k];
        }
        chosenA[(j + 1) * n + j] = 0.0;
        solution[j + 1] -= factor * solution[j];
      }
      // now do the back substitution
      newSolution[n - 1] = solution[n - 1] / chosenA[n * n - 1];
      for (int j = n - 2; j >= 0; j--) {
        double sum = 0.0;
        for (unsigned int k = j + 1; k < n; k++) {
          sum += chosenA[j * n + k] * newSolution[k];
        }
        newSolution[j] = (solution[j] - sum) / chosenA[j * n + j];
      }
      // normalize the solution
      double sum = 0.0;
      for (unsigned int j = 0; j < n; j++) {
        sum += newSolution[j] * newSolution[j];
      }
      sum = sqrt(sum);
      for (unsigned int j = 0; j < n; j++) {
        newSolution[j] /= sum;
      }
      // now copy the new solution to the old solution
      for (unsigned int j = 0; j < n; j++) {
        solution[j] = newSolution[j];
      }
    }
    // now we want to put them back
    for (unsigned int j = 0; j < n; j++) {
      ev[i * n + chosenPositions[j]] = solution[j];
    }
  }
}

template <unsigned int n>
__device__ __forceinline__ void evd(double *A, double *v) {
  static_assert(n==6 || n == 9 || n == 12, "ERROR: n not supported");
  // A is the original matrix
  // after this function, A will be all the eigen vectors
  // v will store all the eigen values
  // first we will do the householder decomposition
  // A = UHU^T
  double U[n * n];
  // set U as identity matrix
  for (unsigned int i = 0; i < n * n; i++) {
    U[i] = 0;
  }
  for (unsigned int i = 0; i < n; i++) {
    U[i * n + i] = 1;
  }
  householderTransformation<n>(
      A, U); // after this step, A_copy will be the tri-diagonal matrix

  // save the tri diagonal system and the U matrix
  double A_tri[n * n];
  for (unsigned int i = 0; i < n * n; i++) {
    A_tri[i] = A[i];
  }

  // now we will do the QR decomposition and jacobi iterations, and only save
  // the eigen values
  qr_tri_diagonal<n>(A, U, false);

  // record the eigen values
  for (unsigned int i = 0; i < n; i++) {
    v[i] = A[i * n + i];
  }

  // now compute the eigen vectors based on the eigen values we got
  computeEigenVectors<n>(A_tri, v, A);

  // we need to check if we got NaN for eigen vectors
  // since the tri diagonal solve does not consider
  // the cases where there's no off diagonal elements
  bool validSolution = true;
  for (unsigned int i = 0; i < n * n; i++) {
    if (A[i] != A[i]) {
      validSolution = false;
      break;
    }
  }

  if (validSolution) {
    // A = U * A^T
    double temp[n * n];
    for (unsigned int i = 0; i < n * n; i++) {
      temp[i] = 0.0;
    }
    for (unsigned int i = 0; i < n; i++) {
      for (unsigned int j = 0; j < n; j++) {
        for (unsigned int k = 0; k < n; k++) {
          temp[i * n + j] += U[i * n + k] * A[j * n + k];
        }
      }
    }
    for (unsigned int i = 0; i < n * n; i++) {
      A[i] = temp[i];
    }
  } else {
    // let's do it the old fashioned way
    qr_tri_diagonal<n>(A_tri, U, true);
    for (unsigned int i = 0; i < n * n; i++) {
      A[i] = U[i];
    }
  }
}
