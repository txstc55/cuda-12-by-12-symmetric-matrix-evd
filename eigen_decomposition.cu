#include <stdio.h>
#include <stdlib.h>
#include <cusolverDn.h>
#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>
#include "householder_n_12.cuh"
// #include "qr_n_12.cuh"
#include "qr_n_12_tri_diagonal.cuh"
#include "evd_12.cuh"
#include <vector>

// Error checking macro
#define CHECK_CUDA(call) { \
    const cudaError_t error = call; \
    if (error != cudaSuccess) { \
        printf("Error: %s:%d, ", __FILE__, __LINE__); \
        printf("code:%d, reason: %s\n", error, cudaGetErrorString(error)); \
        exit(1); \
    } \
}

#define CHECK_CUSOLVER(call) { \
    const cusolverStatus_t error = call; \
    if (error != CUSOLVER_STATUS_SUCCESS) { \
        printf("Error: %s:%d, ", __FILE__, __LINE__); \
        printf("CUSOLVER error code: %d\n", error); \
        exit(1); \
    } \
}

__global__ void generateSymmetricMatrices(double *d_A, int n, unsigned long long seed) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n) return;
    int N = 12;

    // Initialize the random number generator
    curandState_t state;
    curand_init(seed, tid, 0, &state);

    // Generate the upper triangle of the matrix
    for (int row = 0; row < N; ++row) {
        for (int col = row; col < N; ++col) {
            double randomValue = curand_uniform_double(&state) * 2000 - 1000.0;
            d_A[tid * N * N + row * N + col] = randomValue; // Upper triangle
            d_A[tid * N * N + col * N + row] = randomValue; // Mirror to lower triangle
        }
    }

}

__global__ void evd(double* d_A, double* d_e, int n){
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n) return;
    evd_12(d_A + tid * 144, d_e + tid * 12);
}

int main() {
    // ===============================================
    // Setup
    // ===============================================
    cusolverDnHandle_t cusolverH = NULL;
    syevjInfo_t syevj_params = NULL;
    cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_VECTOR; // Compute eigenvalues and eigenvectors
    cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;
    int n = 12; // size of each matrix
    int lda = n;
    int batchSize = 1000000;
    double *d_A = NULL; // Device matrix
    double *d_W = NULL; // Device eigenvalues
    int *d_info = NULL; // info on device
    double *d_work = NULL; // workspace on device
    int lwork = 0; // workspace size
    int info_gpu = 0; // info on host
    cudaEvent_t start, stop;
    float milliseconds = 0;
    unsigned int selected_index = 0;

    // Initialize cuSolver
    CHECK_CUSOLVER(cusolverDnCreate(&cusolverH));
    // Allocate memory on device
    CHECK_CUDA(cudaMalloc((void**)&d_A, sizeof(double) * lda * n * batchSize));
    CHECK_CUDA(cudaMalloc((void**)&d_W, sizeof(double) * n * batchSize));
    CHECK_CUDA(cudaMalloc((void**)&d_info, sizeof(int) * batchSize));
    // Setup the execution configuration
    int threadsPerBlock = 32;
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
    unsigned long long seed = 1234ULL;
    // Launch the kernel to generate random symmetric matrices
    generateSymmetricMatrices<<<blocksPerGrid, threadsPerBlock>>>(d_A, batchSize, seed);
    cudaDeviceSynchronize();


    // ===============================================
    // check and save the matrices
    // ===============================================
    std::vector<double> As;
    As.resize(n * n * batchSize);
    cudaMemcpy(As.data(), d_A, sizeof(double) * n * n * batchSize, cudaMemcpyDeviceToHost);
    printf("Random matrix generated and stored on host and device\n");

    // ===============================================
    // Eigen value decomposition with our code
    // ===============================================
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);
    evd<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_W, batchSize);
    cudaDeviceSynchronize();
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    printf("Time for 1 execution of our eigen value decomposition: %f ms\n", milliseconds);

    // ===============================================
    // check the error for our solution
    // ===============================================
    cudaMemcpy(eigenVectors.data(), d_A, sizeof(double) * n * n * batchSize, cudaMemcpyDeviceToHost);
    cudaMemcpy(eigenValues.data(), d_W, sizeof(double) * n * batchSize, cudaMemcpyDeviceToHost);
    for (int i = 0; i < batchSize; i++){
        double temp[n * n];
        double eigenVector[n * n];
        double eigenValue[n];
        for (int j = 0; j < n * n; j++){
            eigenVector[j] = eigenVectors[i * n * n + j];
        }
        for (int j = 0; j < n; j++){
            eigenValue[j] = eigenValues[i * n + j];
        }

        // now perform U^T * V * U
        for (int j = 0; j < n * n; j++){
            temp[j] = 0;
        }
        for (int j = 0; j < n; j++){
            for (int k = 0; k < n; k++){
                for (int l = 0; l < n; l++){
                    temp[j * n + k] += eigenVector[j * n + l] * eigenValue[l] * eigenVector[k * n + l];
                }
            }
        }
        double diff = 0;
        for (int j = 0; j < n * n; j++){
            diff += abs(temp[j] - As[i * n * n + j]);
        }
        if (diff > maximumError){
            maximumError = diff;
        }
    }
    printf("Maximum error: %lf\n", maximumError);
    maximumError = 0.0;

    // ===============================================
    // setting up Eigen value decomposition with cuSolver
    // ===============================================
    // first we copy back the original matrix
    cudaMemcpy(d_A, As.data(), sizeof(double) * n * n * batchSize, cudaMemcpyHostToDevice);

    // Create syevj parameters
    CHECK_CUSOLVER(cusolverDnCreateSyevjInfo(&syevj_params));

    // Set syevj parameters
    // For example, setting the tolerance and the maximum number of sweeps:
    CHECK_CUSOLVER(cusolverDnXsyevjSetTolerance(syevj_params, 1e-6));
    CHECK_CUSOLVER(cusolverDnXsyevjSetMaxSweeps(syevj_params, 100));

    // Query working space of syevjBatched
    CHECK_CUSOLVER(cusolverDnDsyevjBatched_bufferSize(
        cusolverH,
        jobz,
        uplo,
        n,
        d_A,
        lda,
        d_W,
        &lwork,
        syevj_params,
        batchSize));

    // Allocate workspace for device
    CHECK_CUDA(cudaMalloc((void**)&d_work, sizeof(double) * lwork));
    cudaDeviceSynchronize();


    // ===============================================
    // cusolver execute
    // ===============================================
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);
    // Execute the eigendecomposition
    cusolverDnDsyevjBatched(
        cusolverH,
        jobz,
        uplo,
        n,
        d_A,
        lda,
        d_W,
        d_work,
        lwork,
        d_info,
        syevj_params,
        batchSize);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    printf("Time for 1 execution of cusolver eigen value decomposition: %f ms\n", milliseconds);
    // Synchronize and check for errors
    CHECK_CUDA(cudaDeviceSynchronize());
    CHECK_CUDA(cudaMemcpy(&info_gpu, d_info, sizeof(int), cudaMemcpyDeviceToHost));

    if (info_gpu == 0) {
        printf("Batched syevj execution successful!\n");
    } else {
        printf("Batched syevj execution failed: Info = %d\n", info_gpu);
    }

    // ===============================================
    // check the error for cuSolver solution
    // ===============================================

    cudaMemcpy(eigenVectors.data(), d_A, sizeof(double) * n * n * batchSize, cudaMemcpyDeviceToHost);
    cudaMemcpy(eigenValues.data(), d_W, sizeof(double) * n * batchSize, cudaMemcpyDeviceToHost);

    
    for (int i = 0; i < batchSize; i++){
        double temp[n * n];
        double eigenVector[n * n];
        double eigenValue[n];
        for (int j = 0; j < n * n; j++){
            eigenVector[j] = eigenVectors[i * n * n + j];
        }
        for (int j = 0; j < n; j++){
            eigenValue[j] = eigenValues[i * n + j];
        }

        // now perform U^T * V * U
        for (int j = 0; j < n * n; j++){
            temp[j] = 0;
        }
        for (int j = 0; j < n; j++){
            for (int k = 0; k < n; k++){
                for (int l = 0; l < n; l++){
                    temp[j * n + k] += eigenVector[l * n + j] * eigenValue[l] * eigenVector[l * n + k];
                }
            }
        }
        double diff = 0;
        for (int j = 0; j < n * n; j++){
            diff += abs(temp[j] - As[i * n * n + j]);
        }
        if (diff > maximumError){
            maximumError = diff;
        }
    }
    printf("Maximum error: %lf\n", maximumError);



    // Cleanup
    CHECK_CUDA(cudaFree(d_A));
    CHECK_CUDA(cudaFree(d_W));
    CHECK_CUDA(cudaFree(d_info));
    CHECK_CUDA(cudaFree(d_work));
    CHECK_CUSOLVER(cusolverDnDestroySyevjInfo(syevj_params));
    CHECK_CUSOLVER(cusolverDnDestroy(cusolverH));
    CHECK_CUDA(cudaDeviceReset());

    return EXIT_SUCCESS;
}