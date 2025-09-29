// vmas_cuda.cu
#include <cuda_runtime.h>
#include <cstdio>

__global__ void null_kernel(){ /* no-op */ }

extern "C" void cgogn_cuda_sanity_call() {
    null_kernel<<<1,1>>>();
    cudaError_t e = cudaDeviceSynchronize();
    if (e != cudaSuccess) {
        std::fprintf(stderr, "[CUDA] sanity failed: %s\n", cudaGetErrorString(e));
    } else {
        std::printf("[CUDA] sanity ok\n");
    }
}
