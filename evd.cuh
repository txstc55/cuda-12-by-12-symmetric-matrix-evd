#pragma once
#include "atda_12.cuh"
#include <cuda_runtime.h>
#include <Eigen/Dense>

template <unsigned int N>
__device__ void computeEigenVectorsFromTriDiagonal(const double diagonal[N],
                                                   const double offDiagonal[N],
                                                   const double v[N],
                                                   double *ev) {

  double norm; // For normalization of vectors
  double sum;
  double factor;
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
      bool continueSolve = true;
      for (int k = 0; k < N; k++) {
        if (abs(diagonalCopy[k]) < 1e-16) {
          for (unsigned int j = 0; j < N; j++) {
            ev[i * N + j] = 0;
          }
          ev[i * N + k] = 1;
          continueSolve = false;
          break;
        }
      }
      if (continueSolve) {
        // perform two power iterations
        // we first do the gauss elimination
        // do the first row
        factor =
            offDiagonal[0] / (diagonalCopy[0] == 0 ? 1.0 : diagonalCopy[0]);
        diagonalCopy[1] -= factor * offDiagonal[0];
        solution[1] -= factor * solution[0];
        // do the rest of the rows
        for (int j = 1; j < N - 2; j++) {
          factor =
              offDiagonal[j] / (diagonalCopy[j] == 0 ? 1.0 : diagonalCopy[j]);
          diagonalCopy[j + 1] -= factor * offDiagonal[j];
          solution[j + 1] -= factor * solution[j];
        }

        // do the last row
        factor = offDiagonal[N - 2] /
                 (diagonalCopy[N - 2] == 0 ? 1.0 : diagonalCopy[N - 2]);
        diagonalCopy[N - 1] -= factor * offDiagonal[N - 2];
        solution[N - 1] -= factor * solution[N - 2];

        // now we do the back substitution
        // do the last row
        solution[N - 1] =
            solution[N - 1] /
            (diagonalCopy[N - 1] == 0 ? 1.0 : diagonalCopy[N - 1]);
        // do the rest of the rows
        // #pragma unroll
        for (int j = N - 2; j >= 0; j--) {
          solution[j] = (solution[j] - offDiagonal[j] * solution[j + 1]) /
                        diagonalCopy[j];
        }

        // normalize the result
        sum = 0.0;
        // #pragma unroll
        for (unsigned int j = 0; j < N; j++) {
          sum += solution[j] * solution[j];
        }
        norm = sqrt(sum);
        // #pragma unroll
        for (unsigned int j = 0; j < N; j++) {
          ev[i * N + j] = solution[j] / (norm == 0 ? 1.0 : norm);
        }
      }
    }
  }
}

template <unsigned int N> __device__ void project_to_spd(double *A, double *v) {
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
  atda_12(m_eivalues.data(), MB.data(), A);

  // now we extract the H matrix from tri diagonalization

  for (int i = 0; i < N; i++) {
    v[i] = m_eivalues.data()[i];
  }
}