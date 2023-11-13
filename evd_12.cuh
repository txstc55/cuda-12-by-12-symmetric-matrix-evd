#pragma once
#include <cuda_runtime.h>
#include "householder_n_12.cuh"
#include "qr_n_12_tri_diagonal.cuh"
// this is an attempt to solve for eigen vectors directly from the tri diagonal matrix
// it's working and not working yet
// the problem is the large error it got from the substitution
// __device__ __forceinline__ void computeEigenVectors(const double A[144], const double v[12], double* ev){
//     const int n = 12;
//     double maxAbsVal = 0.0;
//     for (unsigned int i = 0; i < n * n; i++){
//         if (fabs(A[i]) > maxAbsVal){
//             maxAbsVal = fabs(A[i]);
//         }
//     }
//     double scaleFactor = 1.0 / maxAbsVal;

//     for (int i = 0; i < n; i++){
//         double minimumError = 0.0;
//         double lambda = v[i] * scaleFactor;
//         double AForward[n * n];
//         for (unsigned int j = 0; j < n * n; j++){
//             AForward[j] = A[j] * scaleFactor;
//         }
//         for (int j = 0; j < n; j++){
//             AForward[j * n + j] -= lambda;
//         }


//         unsigned int positionsForward[n];
//         for (int j = 0; j < n; j++){
//             positionsForward[j] = j;
//         }
//         // we will implement a row swapping gauss elimination
//         {
//             // deal with the first row
//             if (fabs(AForward[0]) < fabs(AForward[n])){
//                 // swap the first row with the second row
//                 for (int j = 0; j < n; j++){
//                     double temp = AForward[j];
//                     AForward[j] = AForward[j + n];
//                     AForward[j + n] = temp;
//                 }
//                 // int temp = positionsForward[0];
//                 // positionsForward[0] = positionsForward[1];
//                 // positionsForward[1] = temp;
//             }
//             // do the gauss elimination
//             double factor = AForward[n] / AForward[0];
//             for (int j = 0; j < n; j++){
//                 AForward[j + n] -= factor * AForward[j];
//             }
//         }

//         for (int j = 1; j < n - 1; j++){
//             if (fabs(AForward[j * n + j]) < fabs(AForward[(j + 1) * n + j])){
//                 // swap the rows, since it is tri diagonal we are cool with just that
//                 for (int k = 0; k < n; k++){
//                     double temp = AForward[j * n + k];
//                     AForward[j * n + k] = AForward[(j + 1) * n + k];
//                     AForward[(j + 1) * n + k] = temp;
//                 }
//                 // record the swap of positionsForward
//                 // int temp = positionsForward[j];
//                 // positionsForward[j] = positionsForward[j+1];
//                 // positionsForward[j+1] = temp;
//             }
//             // now perform the gauss elimination
//             double factor = AForward[(j + 1) * n + j] / AForward[j * n + j];
//             for (int k = 0; k < n; k++){
//                 AForward[(j + 1) * n + k] -= factor * AForward[j * n + k];
//             }
//         }
//         // printf("Error with forward: %f\n", AForward[n * n - 1]);
//         minimumError = fabs(AForward[n * n - 1]);


//         // now try the backward version
//         double ABackward[n * n];
//         unsigned int positionsBackward[n];
//         for (unsigned int j = 0; j < n * n; j++){
//             ABackward[j] = A[n * n - 1 - j] * scaleFactor;
//         }
//         for (int j = 0; j < n; j++){
//             ABackward[j * n + j] -= lambda;
//         }
//         for (int j = 0; j < n; j++){
//             positionsBackward[j] = n - 1 - j;
//         }
//         // we will implement a row swapping gauss elimination
//         {
//             // deal with the first row
//             if (fabs(ABackward[0]) < fabs(ABackward[n])){
//                 // swap the first row with the second row
//                 for (int j = 0; j < n; j++){
//                     double temp = ABackward[j];
//                     ABackward[j] = ABackward[j + n];
//                     ABackward[j + n] = temp;
//                 }
//                 // int temp = positionsBackward[0];
//                 // positionsBackward[0] = positionsBackward[1];
//                 // positionsBackward[1] = temp;
//             }
//             // do the gauss elimination
//             double factor = ABackward[n] / ABackward[0];
//             for (int j = 0; j < n; j++){
//                 ABackward[j + n] -= factor * ABackward[j];
//             }
//         }

//         for (int j = 1; j < n - 1; j++){
//             if (fabs(ABackward[j * n + j]) < fabs(ABackward[(j + 1) * n + j])){
//                 // swap the rows, since it is tri diagonal we are cool with just that
//                 for (int k = 0; k < n; k++){
//                     double temp = ABackward[j * n + k];
//                     ABackward[j * n + k] = ABackward[(j + 1) * n + k];
//                     ABackward[(j + 1) * n + k] = temp;
//                 }
//                 // record the swap of positionsBackward
//                 // int temp = positionsBackward[j];
//                 // positionsBackward[j] = positionsBackward[j+1];
//                 // positionsBackward[j+1] = temp;
//             }
//             // now perform the gauss elimination
//             double factor = ABackward[(j + 1) * n + j] / ABackward[j * n + j];
//             for (int k = 0; k < n; k++){
//                 ABackward[(j + 1) * n + k] -= factor * ABackward[j * n + k];
//             }
//         }
//         // printf("Error with backward: %f\n", ABackward[n * n - 1]);
//         double chosenA[n * n];
//         unsigned int chosenPositions[n];
//         if (fabs(ABackward[n * n - 1]) < minimumError){
//             minimumError = fabs(ABackward[n * n - 1]);
//             for (unsigned int j = 0; j < n * n; j++){
//                 chosenA[j] = ABackward[j];
//             }
//             for (unsigned int j = 0; j < n; j++){
//                 chosenPositions[j] = positionsBackward[j];
//             }
//         }else{
//             for (unsigned int j = 0; j < n * n; j++){
//                 chosenA[j] = AForward[j];
//             }
//             for (unsigned int j = 0; j < n; j++){
//                 chosenPositions[j] = positionsForward[j];
//             }
//         }

//         // ok actually do the back substitution
//         double solution[n];
//         solution[n - 1] = 1.0;
//         for (int j = n - 2; j >= 0; j--){
//             double sum = 0.0;
//             for (unsigned int k = j + 1; k < n; k++){
//                 sum += chosenA[j * n + k] * solution[k];
//             }
//             solution[j] = -sum / chosenA[j * n + j];
//         }
//         // normalize the solution
//         double sum = 0.0;
//         for (unsigned int j = 0; j < n; j++){
//             sum += solution[j] * solution[j];
//         }
//         sum = sqrt(sum);
//         for (unsigned int j = 0; j < n; j++){
//             solution[j] /= sum;
//         }

//         // now we want to put them back
//         for (unsigned int j = 0; j < n; j++){
//             ev[i * n + chosenPositions[j]] = solution[j];
//         }

//     }
// }


__device__ __forceinline__ void evd_12(double* A, double* v){
    // A is the original matrix
    // after this function, A will be all the eigen vectors
    // v will store all the eigen values
    // first we will do the householder decomposition
    // A = UHU^T
    double U[144];
    // set U as identity matrix
    for (unsigned int i = 0; i < 144; i++){
        U[i] = 0;
    }
    for (unsigned int i = 0; i < 12; i++){
        U[i * 12 + i] = 1;
    }
    householderTransformation(A, U); // after this step, A_copy will be the tri-diagonal matrix
    qr_12_tri_diagonal(A, U);
    for (unsigned int i = 0; i < 12; i++){
        v[i] = A[i * 12 + i];
    }
    for (unsigned int i = 0; i < 144; i++){
        A[i] = U[i];
    }
}
