#include "atda_12.cuh"
#include "evd.cuh"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <ctime>
#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>
#include <cusolverDn.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
using namespace Eigen;

// Error checking macro
#define CHECK_CUDA(call)                                                       \
  {                                                                            \
    const cudaError_t error = call;                                            \
    if (error != cudaSuccess) {                                                \
      printf("Error: %s:%d, ", __FILE__, __LINE__);                            \
      printf("code:%d, reason: %s\n", error, cudaGetErrorString(error));       \
      exit(1);                                                                 \
    }                                                                          \
  }

#define CHECK_CUSOLVER(call)                                                   \
  {                                                                            \
    const cusolverStatus_t error = call;                                       \
    if (error != CUSOLVER_STATUS_SUCCESS) {                                    \
      printf("Error: %s:%d, ", __FILE__, __LINE__);                            \
      printf("CUSOLVER error code: %d\n", error);                              \
      exit(1);                                                                 \
    }                                                                          \
  }

template <unsigned int N>
__global__ void generateSymmetricMatrices(double *d_A, int n,
                                          unsigned long long seed) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= n)
    return;

  // Initialize the random number generator
  curandState_t state;
  curand_init(seed + tid, tid, 0, &state);

  // // Generate the upper triangle of the matrix
  // for (int row = 0; row < N; ++row) {
  //   for (int col = row; col < N; ++col) {
  //     double randomValue = curand_uniform_double(&state) * 2000 - 1000.0;
  //     d_A[tid * N * N + row * N + col] = randomValue; // Upper triangle
  //     d_A[tid * N * N + col * N + row] =
  //         randomValue; // Mirror to lower triangle
  //   }
  // }

  //   // use this code to generate only tri diagonal matrix
  //   for (int row = 0; row < N; ++row) {
  //     d_A[tid * N * N + row * N + row] =
  //         curand_uniform_double(&state) * 2000 - 1000.0;
  //     for (int col = row; col < min(N, row + 1); ++col) {
  //       double randomValue = curand_uniform_double(&state) * 2000 - 1000.0;
  //       d_A[tid * N * N + row * N + col] = randomValue; // Upper triangle
  //       d_A[tid * N * N + col * N + row] =
  //           randomValue; // Mirror to lower triangle
  //     }
  //   }
  // use this code to generate positive definite matrix
  for (int row = 0; row < N; ++row) {
    for (int col = row; col < N; ++col) {
      double randomValue = curand_uniform_double(&state) * 2000 - 1000.0;
      d_A[tid * N * N + row * N + col] = randomValue; // Upper triangle
      d_A[tid * N * N + col * N + row] =
          randomValue; // Mirror to lower triangle
    }
  }
  // perform A = A^T * A
  double temp[N * N];
  for (int i = 0; i < N * N; i++) {
    temp[i] = 0;
  }
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      for (int k = 0; k < N; k++) {
        temp[i * N + j] +=
            d_A[tid * N * N + i * N + k] * d_A[tid * N * N + k * N + j];
      }
    }
  }
  for (int i = 0; i < N * N; i++) {
    d_A[tid * N * N + i] = temp[i];
  }
}

template <typename Scalar, int size>
__device__ __host__ void makePD(Eigen::Matrix<Scalar, size, size> &symMtr,
                                double *d_e) {
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Scalar, size, size>> eigenSolver(
      symMtr);
  for (int i = 0; i < size; i++) {
    d_e[i] = eigenSolver.eigenvalues()[i];
  }
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      symMtr(i, j) = eigenSolver.eigenvectors()(i, j);
    }
  }
}

template <unsigned int N>
__global__ void evd_eigen(double *d_A, double *d_e, int n) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= n)
    return;
  Eigen::Matrix<double, N, N> symMtr;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      symMtr(i, j) = d_A[tid * N * N + i * N + j];
    }
  }

  makePD(symMtr, d_e + tid * N);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      d_A[tid * N * N + i * N + j] = symMtr(i, j);
    }
  }
}

template <unsigned int N>
__global__ void evd_eigen_modified(double *d_A, double *d_e, int n) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= n)
    return;

  project_to_spd<N>(d_A + tid * N * N, d_e + tid * N);
}

int main() {
  // ===============================================
  // Setup
  // ===============================================
  cusolverDnHandle_t cusolverH = NULL;
  syevjInfo_t syevj_params = NULL;
  cusolverEigMode_t jobz =
      CUSOLVER_EIG_MODE_VECTOR; // Compute eigenvalues and eigenvectors
  cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;
  const unsigned int n = 12; // size of each matrix
  int lda = n;
  int batchSize = 100000;
  double *d_A = NULL;    // Device matrix
  double *d_W = NULL;    // Device eigenvalues
  int *d_info = NULL;    // info on device
  double *d_work = NULL; // workspace on device
  int lwork = 0;         // workspace size
  int info_gpu = 0;      // info on host
  cudaEvent_t start, stop;
  float milliseconds = 0;
  unsigned int selected_index = 0;
  std::vector<double> selected_matrix;
  selected_matrix.resize(n * n);
  std::vector<double> selected_eigen_vectors;
  selected_eigen_vectors.resize(n * n);
  std::vector<double> selected_eigen_values;
  selected_eigen_values.resize(n);

  // Initialize cuSolver
  CHECK_CUSOLVER(cusolverDnCreate(&cusolverH));
  // Allocate memory on device
  CHECK_CUDA(cudaMalloc((void **)&d_A, sizeof(double) * lda * n * batchSize));
  CHECK_CUDA(cudaMalloc((void **)&d_W, sizeof(double) * n * batchSize));
  CHECK_CUDA(cudaMalloc((void **)&d_info, sizeof(int) * batchSize));
  // Setup the execution configuration
  int threadsPerBlock = 16;
  int blocksPerGrid = (batchSize + threadsPerBlock - 1) / threadsPerBlock;
  std::vector<double> eigenVectors;
  std::vector<double> eigenValues;
  eigenVectors.resize(n * n * batchSize);
  eigenValues.resize(n * batchSize);
  double maximumError = 0.0;

  // ===============================================
  // Generate random symmetric matrices
  // ===============================================
  // Seed for the random number generator
  // unsigned long long seed = time(NULL);
  unsigned long long seed = 0;
  // Launch the kernel to generate random symmetric matrices
  generateSymmetricMatrices<n>
      <<<blocksPerGrid, threadsPerBlock>>>(d_A, batchSize, seed);
  cudaDeviceSynchronize();

  // ===============================================
  // check and save the matrices
  // ===============================================
  std::vector<double> As;
  As.resize(n * n * batchSize);
  cudaMemcpy(As.data(), d_A, sizeof(double) * n * n * batchSize,
             cudaMemcpyDeviceToHost);
  printf("Random matrix generated and stored on host and device\n");
  cudaMemcpy(selected_matrix.data(), d_A + selected_index * n * n,
             sizeof(double) * n * n, cudaMemcpyDeviceToHost);

  printf("The selected matrix:\n");
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      printf("%lf, ", selected_matrix[i * n + j]);
    }
    printf("\n");
  }
  printf("\n");

  // ===============================================
  // Eigen value decomposition with Eigen code
  // ===============================================
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);
  evd_eigen<n><<<blocksPerGrid, threadsPerBlock>>>(d_A, d_W, batchSize);
  cudaDeviceSynchronize();
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  milliseconds = 0;
  cudaEventElapsedTime(&milliseconds, start, stop);
  printf("Time for 1 execution of our eigen value decomposition: %f ms\n",
         milliseconds);

  // ===============================================
  // check the error for Eigen's solution
  // ===============================================
  cudaMemcpy(eigenVectors.data(), d_A, sizeof(double) * n * n * batchSize,
             cudaMemcpyDeviceToHost);
  cudaMemcpy(eigenValues.data(), d_W, sizeof(double) * n * batchSize,
             cudaMemcpyDeviceToHost);

  for (int i = 0; i < batchSize; i++) {
    double temp[n * n];
    double eigenVector[n * n];
    double eigenValue[n];
    for (int j = 0; j < n * n; j++) {
      eigenVector[j] = eigenVectors[i * n * n + j];
    }
    for (int j = 0; j < n; j++) {
      eigenValue[j] = eigenValues[i * n + j];
    }

    // now perform U^T * V * U
    for (int j = 0; j < n * n; j++) {
      temp[j] = 0;
    }
    for (int j = 0; j < n; j++) {
      for (int k = 0; k < n; k++) {
        for (int l = 0; l < n; l++) {
          temp[j * n + k] +=
              eigenVector[j * n + l] * eigenValue[l] * eigenVector[k * n + l];
          if (temp != temp) {
            printf("NAN detected in eigen vector or eigen values\n");
            exit(1);
          }
        }
      }
    }
    double diff = 0;
    for (int j = 0; j < n * n; j++) {
      diff += abs(temp[j] - As[i * n * n + j]);
    }
    if (diff > maximumError) {
      maximumError = diff;
    }
  }
  printf("Maximum error: %.12lf\n", maximumError);
  maximumError = 0.0;

  cudaMemcpy(selected_eigen_vectors.data(), d_A + selected_index * n * n,
             sizeof(double) * n * n, cudaMemcpyDeviceToHost);
  cudaMemcpy(selected_eigen_values.data(), d_W + selected_index * n,
             sizeof(double) * n, cudaMemcpyDeviceToHost);
  printf("The selected matrix's eigen vectors using Eigen's method:\n");
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      printf("%lf, ", selected_eigen_vectors[i * n + j]);
    }
    printf("\n");
  }
  printf("\n");
  printf("The selected matrix's eigen values using Eigen's method:\n");
  for (int i = 0; i < n; i++) {
    printf("%lf, ", selected_eigen_values[i]);
  }
  printf("\n\n");

  // ===============================================
  // Eigen value decomposition with our code
  // ===============================================
  // first we copy back the original matrix
  cudaMemcpy(d_A, As.data(), sizeof(double) * n * n * batchSize,
             cudaMemcpyHostToDevice);
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);
  evd_eigen_modified<n>
      <<<blocksPerGrid, threadsPerBlock>>>(d_A, d_W, batchSize);
  cudaDeviceSynchronize();
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  milliseconds = 0;
  cudaEventElapsedTime(&milliseconds, start, stop);
  printf("Time for 1 execution of our eigen value decomposition: %f ms\n",
         milliseconds);
  // ===============================================
  // check the error for our solution
  // ===============================================
  cudaMemcpy(eigenVectors.data(), d_A, sizeof(double) * n * n * batchSize,
             cudaMemcpyDeviceToHost);
  cudaMemcpy(eigenValues.data(), d_W, sizeof(double) * n * batchSize,
             cudaMemcpyDeviceToHost);

  for (int i = 0; i < batchSize; i++) {
    double temp[n * n];
    double eigenVector[n * n];
    double eigenValue[n];
    for (int j = 0; j < n * n; j++) {
      eigenVector[j] = eigenVectors[i * n * n + j];
    }
    for (int j = 0; j < n; j++) {
      eigenValue[j] = eigenValues[i * n + j];
    }

    // now perform U^T * V * U
    for (int j = 0; j < n * n; j++) {
      temp[j] = 0;
    }
    for (int j = 0; j < n; j++) {
      for (int k = 0; k < n; k++) {
        for (int l = 0; l < n; l++) {
          temp[j * n + k] +=
              eigenVector[j * n + l] * eigenValue[l] * eigenVector[k * n + l];
          if (temp != temp) {
            printf("NAN detected in eigen vector or eigen values\n");
            exit(1);
          }
        }
      }
    }
    double diff = 0;
    for (int j = 0; j < n * n; j++) {
      diff += abs(temp[j] - As[i * n * n + j]);
    }
    if (diff > maximumError) {
      maximumError = diff;
    }
  }
  printf("Maximum error: %.12lf\n", maximumError);
  maximumError = 0.0;

  cudaMemcpy(selected_eigen_vectors.data(), d_A + selected_index * n * n,
             sizeof(double) * n * n, cudaMemcpyDeviceToHost);
  cudaMemcpy(selected_eigen_values.data(), d_W + selected_index * n,
             sizeof(double) * n, cudaMemcpyDeviceToHost);
  printf("The selected matrix's eigen vectors using our method:\n");
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      printf("%lf, ", selected_eigen_vectors[i * n + j]);
    }
    printf("\n");
  }
  printf("\n");
  printf("The selected matrix's eigen values using our method:\n");
  for (int i = 0; i < n; i++) {
    printf("%lf, ", selected_eigen_values[i]);
  }
  printf("\n\n");

  // ===============================================
  // setting up Eigen value decomposition with cuSolver
  // ===============================================
  // first we copy back the original matrix
  cudaMemcpy(d_A, As.data(), sizeof(double) * n * n * batchSize,
             cudaMemcpyHostToDevice);

  // Create syevj parameters
  CHECK_CUSOLVER(cusolverDnCreateSyevjInfo(&syevj_params));

  // Set syevj parameters
  // For example, setting the tolerance and the maximum number of sweeps:
  CHECK_CUSOLVER(cusolverDnXsyevjSetTolerance(syevj_params, 1e-6));
  CHECK_CUSOLVER(cusolverDnXsyevjSetMaxSweeps(syevj_params, 1000));

  // Query working space of syevjBatched
  CHECK_CUSOLVER(cusolverDnDsyevjBatched_bufferSize(cusolverH, jobz, uplo, n,
                                                    d_A, lda, d_W, &lwork,
                                                    syevj_params, batchSize));

  // Allocate workspace for device
  CHECK_CUDA(cudaMalloc((void **)&d_work, sizeof(double) * lwork));
  cudaDeviceSynchronize();

  // ===============================================
  // cusolver execute
  // ===============================================
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);
  // Execute the eigendecomposition
  cusolverDnDsyevjBatched(cusolverH, jobz, uplo, n, d_A, lda, d_W, d_work,
                          lwork, d_info, syevj_params, batchSize);
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  milliseconds = 0;
  cudaEventElapsedTime(&milliseconds, start, stop);
  printf("Time for 1 execution of cusolver eigen value decomposition: %f ms\n",
         milliseconds);
  // Synchronize and check for errors
  CHECK_CUDA(cudaDeviceSynchronize());
  CHECK_CUDA(
      cudaMemcpy(&info_gpu, d_info, sizeof(int), cudaMemcpyDeviceToHost));

  if (info_gpu == 0) {
    printf("Batched syevj execution successful!\n");
  } else {
    printf("Batched syevj execution failed: Info = %d\n", info_gpu);
  }

  // ===============================================
  // check the error for cuSolver solution
  // ===============================================

  cudaMemcpy(eigenVectors.data(), d_A, sizeof(double) * n * n * batchSize,
             cudaMemcpyDeviceToHost);
  cudaMemcpy(eigenValues.data(), d_W, sizeof(double) * n * batchSize,
             cudaMemcpyDeviceToHost);

  maximumError = 0.0;
  for (int i = 0; i < batchSize; i++) {
    double temp[n * n];
    double eigenVector[n * n];
    double eigenValue[n];
    for (int j = 0; j < n * n; j++) {
      eigenVector[j] = eigenVectors[i * n * n + j];
    }
    for (int j = 0; j < n; j++) {
      eigenValue[j] = eigenValues[i * n + j];
    }

    // now perform U^T * V * U
    for (int j = 0; j < n * n; j++) {
      temp[j] = 0;
    }
    for (int j = 0; j < n; j++) {
      for (int k = 0; k < n; k++) {
        for (int l = 0; l < n; l++) {
          temp[j * n + k] +=
              eigenVector[l * n + j] * eigenValue[l] * eigenVector[l * n + k];
        }
      }
    }
    double diff = 0;
    for (int j = 0; j < n * n; j++) {
      diff += abs(temp[j] - As[i * n * n + j]);
    }
    if (diff > maximumError) {
      maximumError = diff;
    }
  }
  printf("Maximum error: %.12lf\n", maximumError);
  cudaMemcpy(selected_eigen_vectors.data(), d_A + selected_index * n * n,
             sizeof(double) * n * n, cudaMemcpyDeviceToHost);
  cudaMemcpy(selected_eigen_values.data(), d_W + selected_index * n,
             sizeof(double) * n, cudaMemcpyDeviceToHost);
  printf("The selected matrix's eigen vectors using cusolver:\n");
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      printf("%lf, ", selected_eigen_vectors[j * n + i]);
    }
    printf("\n");
  }
  printf("\n");
  printf("The selected matrix's eigen values using cusolver:\n");
  for (int i = 0; i < n; i++) {
    printf("%lf, ", selected_eigen_values[i]);
  }
  printf("\n\n");

  // Cleanup
  CHECK_CUDA(cudaFree(d_A));
  CHECK_CUDA(cudaFree(d_W));
  CHECK_CUDA(cudaFree(d_info));
  CHECK_CUDA(cudaFree(d_work));
  CHECK_CUDA(cudaEventDestroy(start));
  CHECK_CUDA(cudaEventDestroy(stop));
  CHECK_CUSOLVER(cusolverDnDestroySyevjInfo(syevj_params));
  CHECK_CUSOLVER(cusolverDnDestroy(cusolverH));
  CHECK_CUDA(cudaDeviceReset());

  return EXIT_SUCCESS;
}