// vmas_cuda.cu
#include "vmas_cuda_shared.h"
#include <cgogn/geometry/types/cuda_plain_type.h>
#include <cstdio>
#include <cuda_runtime.h>

using PlainSphericalQuadric = cgogn::cuda::PlainSphericalQuadric;

__device__ inline float3 load_float3(const float* arr, int idx)
{
	return make_float3(arr[3 * idx], arr[3 * idx + 1], arr[3 * idx + 2]);
}
__device__ inline float4 load_float4(const float* arr, int idx)
{
	return make_float4(arr[4 * idx], arr[4 * idx + 1], arr[4 * idx + 2], arr[4 * idx + 3]);
}

__device__ inline float dot_float3(const float3& a, const float3& b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

__device__ inline float compute_distance_term(const float3& samplePos, const float3& projPos,
											  const float3& sphereCenter, float sphereRadius, int distanceMode)
{
	switch (distanceMode)
	{
	case 0: { // SPHERE_EUCLIDEAN_DISTANCE
		const float3 diff =
			make_float3(projPos.x - sphereCenter.x, projPos.y - sphereCenter.y, projPos.z - sphereCenter.z);
		const float d = sqrtf(dot_float3(diff, diff)) - sphereRadius;
		return d * d;
	}
	case 1: { // SPHERE_CENTER_DISTANCE
		const float3 diff =
			make_float3(samplePos.x - sphereCenter.x, samplePos.y - sphereCenter.y, samplePos.z - sphereCenter.z);
		const float d = sqrtf(dot_float3(diff, diff));
		return d * d;
	}
	case 2: { // SPHERE_POWER_DISTANCE
		const float3 diff =
			make_float3(samplePos.x - sphereCenter.x, samplePos.y - sphereCenter.y, samplePos.z - sphereCenter.z);
		return dot_float3(diff, diff) - sphereRadius * sphereRadius;
	}
	case 3: { // SPHERE_POWER_DISTANCE_SQUARED
		const float3 diff =
			make_float3(samplePos.x - sphereCenter.x, samplePos.y - sphereCenter.y, samplePos.z - sphereCenter.z);
		const float pd = dot_float3(diff, diff) - sphereRadius * sphereRadius;
		return pd * pd;
	}
	default:
		return 0.f;
	}
}

__device__ inline float compute_SQEM_energy(const float3& samplePos, const float3& projPos, const float3& proNormal,
											const float* quadricA, const float* quadricB, const float quadricC,
											const float3& sphereCenter, float sphereRadius)
{
	float vec[4] = {sphereCenter.x, sphereCenter.y, sphereCenter.z, sphereRadius};

	float sqem = 0.f;
	for (int r = 0; r < 4; ++r)
	{
		float Av = 0.f;
		for (int c = 0; c < 4; ++c)
			Av += quadricA[r * 4 + c] * vec[c];
		sqem += 0.5f * vec[r] * Av;
	}
	float dot_b = 0.f;
	for (int i = 0; i < 4; i++)
	{
		dot_b += quadricB[i] * vec[i];
	}
	sqem -= dot_b;
	sqem += quadricC;
	return sqem;
}

__device__ inline float compute_energy(const float3& samplePos, const float3& projPos, const float3& proNormal,
									   const PlainSphericalQuadric& quadric, const float3& sphereCenter,
									   float sphereRadius, float lambda, int distanceMode)
{
	const float distanceTerm = compute_distance_term(samplePos, projPos, sphereCenter, sphereRadius, distanceMode);
	const float sqemTerm =
		compute_SQEM_energy(samplePos, projPos, proNormal, quadric.A, quadric.b, quadric.c, sphereCenter, sphereRadius);
	return lambda * distanceTerm + sqemTerm;
}

__global__ void compute_membership_kernel(
	const int nbSamples, const float4* __restrict__ d_samplesPos, const float4* __restrict__ d_samplesProj,
	const float4* __restrict__ d_samplesNorm, const PlainSphericalQuadric* __restrict__ d_quadric, const int nbSpheres,
	const float4* __restrict__ d_spheresPos, const float* __restrict__ d_spheresRadius, const float lambda,
	const float tau, const float eps, const float thetaGap, const int distanceMode,
	MembershipEntry* __restrict__ d_outEntries, int* __restrict__ d_outCounts)
{

	const int tid = blockIdx.x * blockDim.x + threadIdx.x;
	const int stride = blockDim.x * gridDim.x;

	for (int sampleId = tid; sampleId < nbSamples; sampleId += stride)
	{

		const float3 samplePos = make_float3(d_samplesPos[sampleId].x, d_samplesPos[sampleId].y, d_samplesPos[sampleId].z);
		const float3 sampleProj = make_float3(d_samplesProj[sampleId].x, d_samplesProj[sampleId].y, d_samplesProj[sampleId].z);
		const float3 sampleNorm = make_float3(d_samplesNorm[sampleId].x, d_samplesNorm[sampleId].y, d_samplesNorm[sampleId].z);
		const PlainSphericalQuadric& quadric = d_quadric[sampleId];
		float bestEnergys[VMAS_MAX_TOPK];
		int bestIndices[VMAS_MAX_TOPK];
		int count = 0;
		for (int s = 0; s < nbSpheres; s++)
		{
			const float3 center = make_float3(d_spheresPos[s].x, d_spheresPos[s].y, d_spheresPos[s].z);
			const float radius = d_spheresRadius[s];
			const float energy =
				compute_energy(samplePos, sampleProj, sampleNorm, quadric, center, radius, lambda, distanceMode);
			if (count < VMAS_MAX_TOPK)
			{
				int j = count - 1;
				while (j >= 0 && bestEnergys[j] > energy)
				{
					bestEnergys[j + 1] = bestEnergys[j];
					bestIndices[j + 1] = bestIndices[j];
					j--;
				}
				bestEnergys[j + 1] = energy;
				bestIndices[j + 1] = s;
				++count;
			}
			else if (energy < bestEnergys[count - 1])
			{
				int j = count - 2;
				while (j >= 0 && bestEnergys[j] > energy)
				{
					bestEnergys[j + 1] = bestEnergys[j];
					bestIndices[j + 1] = bestIndices[j];
					j--;
				}
				bestEnergys[j + 1] = energy;
				bestIndices[j + 1] = s;
			}
		}
		if (count == 0)
		{
			d_outCounts[sampleId] = 0;
			continue;
		}

		float logE[VMAS_MAX_TOPK];
		for (int i = 0; i < count; i++)
		{
			logE[i] = logf(bestEnergys[i] + eps);
		}
		float gap[VMAS_MAX_TOPK];
		for (int i = 0; i + 1 < count; ++i)
		{
			gap[i] = logE[i + 1] - logE[i];
		}
		int pike = count - 1;
		for (int i = 0; i + 1 < count; ++i)
		{
			if (gap[i] > thetaGap)
			{
				pike = i;
				break;
			}
		}
		const int keepCount = pike + 1;
		MembershipEntry* dst = &d_outEntries[sampleId * VMAS_MAX_TOPK];

		float Emin = bestEnergys[0];
		float denom = 0.f;
		for (int i = 0; i < keepCount; ++i)
		{
			const float w = expf((-bestEnergys[i]-Emin)/ tau);
			denom += w;
			dst[i].sphere_idx = bestIndices[i];
			dst[i].weight = w;
		}
		if (keepCount == 0)
		{
			dst[0].sphere_idx = bestIndices[0];
			dst[0].weight = 1.f;
			denom = 1.f;
			d_outCounts[sampleId] = 1;
			continue;
		}
		if (denom <= 0.f)
			denom = 1.f;
		for (int i = 0; i < keepCount; ++i)
		{
			dst[i].weight /= denom;
		}
		d_outCounts[sampleId] = keepCount;
	}
}

extern "C" cudaError_t cgogn_compute_membership(const int nbSamples, const float4* d_samplesPos,
												const float4* d_samplesProj, const float4* d_samplesNorm,
												const PlainSphericalQuadric* d_quadric, const int nbSpheres,
												const float4* d_spheresPos, const float* d_spheresRadius,
												const float lambda, const float tau, const float eps,
												const float thetaGap, const int distanceMode,
												MembershipEntry* d_outEntries, int* d_outCounts)
{
	if (nbSamples == 0 || nbSpheres == 0)
		return cudaSuccess;

	cudaDeviceProp prop{};
	cudaGetDeviceProperties(&prop, 0);
	const int sms = prop.multiProcessorCount;
	const int blockSize = 64;
	const int gridSize = sms * 16;

	compute_membership_kernel<<<gridSize, blockSize>>>(nbSamples, d_samplesPos, d_samplesProj, d_samplesNorm, d_quadric,
													   nbSpheres, d_spheresPos, d_spheresRadius, lambda, tau, eps,
													   thetaGap, distanceMode, d_outEntries, d_outCounts);

	return cudaGetLastError();
}