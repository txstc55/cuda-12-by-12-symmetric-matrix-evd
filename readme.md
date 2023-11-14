# 12 By 12 Symmetric Matrix Eigen Value Decomposition

This repo contains code for performing 12 by 12 symmetric matrix's eigen value decomposition.

To use the code, just 

``` c++
include "evd_12.cuh"
__global__ void evd(double* d_A, double* d_e, int n){
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n) return;
    evd_12(d_A + tid * 144, d_e + tid * 12);
}

int main(){
    // ...
    evd<<<BLOCKS, THREADS>>>(d_A, d_e, n);
    // ...
}
```

The results will be the Eigen vectors, which are stored in each column in `d_A`, the eigen values stored in `d_e`. To reconstruct the matrix, just do `d_A * d_e * d_A.T`.

The performance tested on RTX4090 for 1M matrix is 330ms, comparing to CUSOLVER's 580ms, it's almost 1.75x faster. You can also just run the `eigen_decomposition.cu` to see the performance. To compile, run:

``` bash
nvcc eigen_decomposition.cu -o eigen.out -O3 -gencode arch=compute_86,code=sm_86 -use_fast_math -lcusolver
./eigen.out
```

The reason for `arch=compute_86,code=sm_86` is simply because I'm using a 4090, you can replace that part with whatever arch you want. But do notice that the support for double is better for the later archs.

Note that for the result that cusolver gives, the eigen vectors are stored in each row, which is different from what I have here (stored in each column). So you need to transpose the result from cusolver to compare with my result.

## Verification

The binary will print out one of the matrix and its eigen vectors and eigen values produced both by us and cusolver. The error term is computed as the sum of l1 norm of the reconstructed matrix and the original matrix. You can also use the `verify.py` to verify the correctness of the code where you need to plug in the results.

There is also a NaN value check only for my result, so if anything happens it will just quit. 

## Code cannot run

Try to decrease the batchSize and run again. In reality, although 1M matrices should only take around 1GB memory, CUSOLVER actually takes more memory than needed.
