#pragma once
#include <cuda_runtime.h>
#include <cgogn/geometry/types/cuda_plain_type.h>

constexpr int VMAS_MAX_TOPK = 8;
using uint32 = std::uint32_t;
#define CUDA_OK(expr, msg)                                                                                     \
    do                                                                                                         \
    {                                                                                                          \
        cudaError_t _e = (expr);                                                                               \
        if (_e != cudaSuccess)                                                                                 \
        {                                                                                                      \
            std::fprintf(stderr, "CUDA error (%s): %s at %s:%d\n",                                             \
                        (msg), cudaGetErrorString(_e), __FILE__, __LINE__);                                     \
            std::exit(EXIT_FAILURE);                                                                           \
        }                                                                                                      \
    } while (0)
struct MembershipEntry
{
	uint32 sphere_idx;
	float weight;
};

extern "C" cudaError_t cgogn_compute_membership(int nbSamples, const float4* d_samplesPos, const float4* d_samplesProj,
												const float4* d_samplesNorm, const cgogn::cuda::PlainSphericalQuadric* d_quadric, int nbSpheres,
												const float4* d_spheresPos, const float* d_spheresRadius, float lambda,
												float tau, float eps, float thetaGap, int distanceMode,
												MembershipEntry* d_outEntries, int* d_outCounts);
