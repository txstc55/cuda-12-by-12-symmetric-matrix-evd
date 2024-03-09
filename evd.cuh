#pragma once
#include <Eigen/Dense>
#include <cuda_runtime.h>

template <unsigned int N>
__device__ void computeEigenVectorsFromTriDiagonal(const double diagonal[N],
                                                   const double offDiagonal[N],
                                                   const double v[N],
                                                   double *ev) {

  double norm; // For normalization of vectors
  double sum;
  double factor;
  for (int i = 0; i < N; i++) {
    if (abs(diagonal[i] - v[i]) < 1e-16) {
      ev[0] = 0.0 / 0.0; // force NaN
      return;
    }
  }
  for (unsigned int i = 0; i < N; i++) {
    double solution[N] = {0}; // Temporary vector for eigenvector computations
    double lambda = v[i];
    double diagonalCopy[N];
    solution[0] = 1.0;

    for (int iteration = 0; iteration < 2; iteration++) {
      for (int j = 0; j < N; j++) {
        diagonalCopy[j] = diagonal[j] - lambda;
      }
      // do a check
      // perform two power iterations
      // we first do the gauss elimination
      // do the first row
      factor = offDiagonal[0] / diagonalCopy[0];
      diagonalCopy[1] -= factor * offDiagonal[0];
      solution[1] -= factor * solution[0];
      // do the rest of the rows
      for (int j = 1; j < N - 2; j++) {
        factor = offDiagonal[j] / diagonalCopy[j];
        diagonalCopy[j + 1] -= factor * offDiagonal[j];
        solution[j + 1] -= factor * solution[j];
      }

      // do the last row
      factor = offDiagonal[N - 2] / diagonalCopy[N - 2];
      diagonalCopy[N - 1] -= factor * offDiagonal[N - 2];
      solution[N - 1] -= factor * solution[N - 2];

      // now we do the back substitution
      // do the last row
      solution[N - 1] = solution[N - 1] / diagonalCopy[N - 1];
      // do the rest of the rows
      // #pragma unroll
      for (int j = N - 2; j >= 0; j--) {
        solution[j] =
            (solution[j] - offDiagonal[j] * solution[j + 1]) / diagonalCopy[j];
      }

      // normalize the result
      sum = 0.0;
      // #pragma unroll
      for (unsigned int j = 0; j < N; j++) {
        sum += solution[j] * solution[j];
      }
      norm = sqrt(sum);
      if (norm != norm) {
        ev[0] = 0.0 / 0.0; // force NaN
        return;
      }
      // #pragma unroll
      for (unsigned int j = 0; j < N; j++) {
        ev[i * N + j] = solution[j] / norm;
      }
    }
  }
}

__device__ __forceinline__ double dot_product(const double *v1,
                                              const double *v2, int size) {
  double dot = 0.0;
  for (int i = 0; i < size; ++i) {
    dot += v1[i] * v2[i];
  }
  return dot;
}

template <unsigned int n>
__device__ __forceinline__ void householderTransformation(double *A,
                                                          double *U) {
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
__device__ __forceinline__ void find_pivot_and_rotate(double *A, double *E,
                                                      bool updateEigenVectors) {
  const double epsilon =
      1e-12; // A small threshold to handle very small numbers
  for (unsigned int iteration = 0; iteration < 5; iteration++) {
    double max_value =
        epsilon; // Initialize with epsilon to handle small numbers
    int p = -1;
    int q = -1;

    // Find the indices of the largest off-diagonal element in A
    for (int i = 0; i < n; ++i) {
      for (int j = i + 1; j < n; ++j) {
        double abs_value = fabs(A[i * n + j]);
        if (abs_value > max_value) {
          max_value = abs_value;
          p = i;
          q = j;
        }
      }
    }

    if (max_value <= epsilon) { // Check for convergence
      break;
    }

    // Calculate the rotation angle
    double a_pp = A[p * n + p];
    double a_qq = A[q * n + q];
    double a_pq = A[p * n + q];
    double theta =
        0.5 * atan2(2.0 * a_pq,
                    a_qq - a_pp); // Use atan2 for better numerical stability
    double c = cos(theta);
    double s = sin(theta);

    // Perform the rotation
    for (int i = 0; i < n; ++i) {
      if (i != p && i != q) {
        double a_ip = A[i * n + p];
        double a_iq = A[i * n + q];
        A[i * n + p] = c * a_ip - s * a_iq;
        A[i * n + q] = s * a_ip + c * a_iq;
        A[p * n + i] = A[i * n + p]; // Maintain symmetry
        A[q * n + i] = A[i * n + q]; // Maintain symmetry
      }
    }

    // Update the diagonal elements
    A[p * n + p] = c * c * a_pp + s * s * a_qq - 2.0 * s * c * a_pq;
    A[q * n + q] = s * s * a_pp + c * c * a_qq + 2.0 * s * c * a_pq;

    // Zero out the pivot element
    A[p * n + q] = 0.0;
    A[q * n + p] = 0.0; // Maintain symmetry

    if (updateEigenVectors) {
      // Update the rotation matrix
      for (int i = 0; i < n; ++i) {
        double e_ip = E[i * n + p];
        double e_iq = E[i * n + q];
        E[i * n + p] = c * e_ip - s * e_iq;
        E[i * n + q] = s * e_ip + c * e_iq;
      }
    }
  }
}

template <unsigned int n>
__device__ __forceinline__ void
qr_tri_diagonal(double A[n * n], double E[n * n], bool updateEigenVectors) {
  double ACopy[n * n];
  for (int i = 0; i < n * n; ++i) {
    ACopy[i] = A[i];
  }

  double eq[n * n];
  double e_tmp, q_tmp;

  // now we do the jacobi iterations
  bool keepDoingJacobian = true;
  while (keepDoingJacobian) {
    find_pivot_and_rotate<n>(ACopy, E, updateEigenVectors);
    keepDoingJacobian = false;
    for (int i = 0; i < n; i++) {
      for (int j = i + 1; j < n; j++) {
        double absv = fabs(ACopy[i * n + j]);
        if (absv > 1e-6) {
          keepDoingJacobian = true;
          break;
        }
      }
    }
  }

  for (unsigned int i = 0; i < n; i++) {
    A[i * n + i] = ACopy[i * n + i];
  }
}

template <unsigned int n>
__device__ __forceinline__ void evd_ours_small(double *A, double *v) {
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
  double diagonal[n];
  double offDiagonal[n];
  computeEigenVectorsFromTriDiagonal<n>(diagonal, offDiagonal, v, A);

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
    for (unsigned int i = 0; i < n; i++) {
      for (unsigned int j = 0; j < n; j++) {
        A[i * n + j] = temp[j * n + i];
      }
    }
  } else {
    // let's do it the old fashioned way
    qr_tri_diagonal<n>(A_tri, U, true);
    for (unsigned int i = 0; i < n; i++) {
      for (unsigned int j = 0; j < n; j++) {
        A[i * n + j] = U[j * n + i];
      }
    }
  }
}

template <unsigned int N> __device__ void evd_ours(double *A, double *v) {
  // case for N == 3
  if (N == 3) {
    evd_ours_small<N>(A, v);
    return;
  }
  // cases for N <= 6
  if (N <= 6) {
    Eigen::Matrix<double, N, N> symMtr;
    for (int i = 0; i < N * N; i++) {
      symMtr.data()[i] = A[i];
    }
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, N, N>> eigenSolver(
        symMtr);
    Eigen::Matrix<double, N, N> B = eigenSolver.eigenvectors();
    for (int i = 0; i < N; i++) {
      v[i] = eigenSolver.eigenvalues()[i];
    }
    for (int i = 0; i < N * N; i++) {
      A[i] = B.data()[i];
    }
    return;
  }
  // first declare some typedefs
  typedef typename Eigen::internal::plain_col_type<Eigen::Matrix<double, N, N>,
                                                   double>::type VectorType;
  typedef typename Eigen::internal::plain_col_type<Eigen::Matrix<double, N, N>,
                                                   double>::type RealVectorType;
  typename Eigen::Tridiagonalization<
      Eigen::Matrix<double, N, N>>::SubDiagonalType m_subdiag;
  typename Eigen::Tridiagonalization<
      Eigen::Matrix<double, N, N>>::CoeffVectorType m_hcoeffs;
  typedef Eigen::Matrix<double, N, N, Eigen::ColMajor, N, N> EigenvectorsType;
  typedef typename Eigen::Tridiagonalization<Eigen::Matrix<double, N, N>>::
      HouseholderSequenceType HouseholderSequenceType;

  // first we copy matrix A
  Eigen::Matrix<double, N, N> matrix;
  for (int i = 0; i < N * N; i++) {
    matrix.data()[i] = A[i];
  }

  // pre allocate some spaces
  VectorType m_workspace;
  RealVectorType m_eivalues;
  EigenvectorsType m_eivec;

  // declare some aliases
  RealVectorType &diag = m_eivalues;
  EigenvectorsType &mat = m_eivec;

  // rescale the matrix
  mat = matrix.template triangularView<Eigen::Lower>();
  double scale = mat.cwiseAbs().maxCoeff();
  if (scale == 0)
    scale = double(1);
  mat.template triangularView<Eigen::Lower>() /= scale;

  // space allocation
  m_eivalues.resize(N, 1);
  m_workspace.resize(N);
  m_subdiag.resize(N - 1);
  m_hcoeffs.resize(N - 1);
  m_eivalues.resize(N, 1);

  // first we perform the tri diagonalization
  Eigen::internal::tridiagonalization_inplace(mat, m_hcoeffs);
  mat = mat.real();
  diag = mat.diagonal();
  m_subdiag = mat.template diagonal<-1>();

  // copy the diagonal and tri diagonal data
  double diag_data[N];
  for (int i = 0; i < N; i++) {
    diag_data[i] = mat.data()[i * N + i];
  }
  double sub_diag[N - 1];
  for (int i = 0; i < N - 1; i++) {
    sub_diag[i] = mat.data()[(i)*N + i + 1];
  }

  // now we compute the diagonalization, without computing the eigen values
  Eigen::internal::computeFromTridiagonal_impl(m_eivalues, m_subdiag, N, false,
                                               m_eivec);

  // now we check if the eigen values are >=0
  // because eigen values are already sorted
  // we can just check the first one

  if (m_eivalues.data()[0] >= 0) {
    return;
  }

  // now we compute the eigen vectors
  // we perform a tri diagonal solve
  Eigen::Matrix<double, N, N> B;
  computeEigenVectorsFromTriDiagonal<N>(diag_data, sub_diag, m_eivalues.data(),
                                        B.data());
  // tri diagonal solve failed
  if (B.data()[0] != B.data()[0]) {
    Eigen::Matrix<double, N, N> symMtr;
    for (int i = 0; i < N * N; i++) {
      symMtr.data()[i] = A[i];
    }
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, N, N>> eigenSolver(
        symMtr);
    B = eigenSolver.eigenvectors();
    for (int i = 0; i < N; i++) {
      v[i] = eigenSolver.eigenvalues()[i];
    }
    for (int i = 0; i < N * N; i++) {
      A[i] = B.data()[i];
    }
    return;
  }

  // then we extract the householder matrix
  auto house = HouseholderSequenceType(mat, m_hcoeffs)
                   .setLength(mat.rows() - 1)
                   .setShift(1);

  mat.diagonal().setOnes();
  for (int k = N - 2; k >= 0; --k) {
    mat.bottomRightCorner(N - k - 1, N - k - 1)
        .applyHouseholderOnTheLeft(house.essentialVector(k), m_hcoeffs.coeff(k),
                                   m_workspace.data());

    // clear the off diagonal vector
    mat.col(k).tail(N - k - 1).setZero();
  }

  // now we get the eigen values
  Eigen::Matrix<double, N, N> MB = mat * B; // multiply the householder matrix
  for (int i = 0; i < N; i++) {
    m_eivalues.data()[i] *= scale;
  }

  // we put it back
  // atda_12(m_eivalues.data(), MB.data(), A);
  for (int i = 0; i < N * N; i++) {
    A[i] = MB.data()[i];
  }

  // now we extract the H matrix from tri diagonalization

  for (int i = 0; i < N; i++) {
    v[i] = m_eivalues.data()[i];
  }
}