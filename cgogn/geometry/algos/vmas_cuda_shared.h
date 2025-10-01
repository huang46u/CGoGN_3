#pragma once
#include <cuda_runtime.h>
constexpr int VMAS_MAX_TOPK = 8;
#define CUDA_OK(expr)                                                                                                  \
	do                                                                                                                 \
	{                                                                                                                  \
		cudaError_t _e = (expr);                                                                                       \
		if (_e != cudaSuccess)                                                                                         \
		{                                                                                                              \
			std::fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(_e), __FILE__, __LINE__);              \
			std::exit(1);                                                                                              \
		}                                                                                                              \
	} while (0)

struct MembershipEntry
{
	int sphere_idx;
	float weight;
};

extern "C" cudaError_t cgogn_compute_membership(int nbSamples, const float* d_samplesPos, const float* d_samplesProj,
												const float* d_samplesNorm, const float* d_quadricA,
												const float* d_quadricB, const float* d_quadricC, int nbSpheres,
												const float* d_spheresPos, const float* d_spheresRadius, float lambda,
												float tau, float eps, float thetaGap, int distanceMode,
												MembershipEntry* d_outEntries, int* d_outCounts);
