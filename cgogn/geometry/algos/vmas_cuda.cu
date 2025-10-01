// vmas_cuda.cu
#include "vmas_cuda_shared.h"
#include <cstdio>
#include <cuda_runtime.h>

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

__device__ inline float compute_distance_term(const float3& samplePos, const float3& projPos, const float3& sphereCenter,
											  float sphereRadius, int distanceMode)
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
	float Av[4] = {0.f, 0.f, 0.f, 0.f};
	// matrix-vector multiplication Av = A * v
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++)
		{
			Av[r] += quadricA[r * 4 + c] * vec[c];
		}
	}
	float sqem = 0.f;
	for (int i = 0; i < 4; i++)
	{
		sqem += 0.5f * vec[i] * Av[i];
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
									   const float* quadricA, const float* quadricB, const float quadricC,
									   const float3& sphereCenter, float sphereRadius, float lambda, int distanceMode)
{
	const float distanceTerm = compute_distance_term(samplePos, projPos, sphereCenter, sphereRadius, distanceMode);
	const float sqemTerm =
		compute_SQEM_energy(samplePos, projPos, proNormal, quadricA, quadricB, quadricC, sphereCenter, sphereRadius);
	return lambda * distanceTerm + sqemTerm;
}

__global__ void compute_membership_kernel(const int nbSamples,
	const float* d_samplesPos,
	const float* d_samplesProj,
	const float* d_samplesNorm,
	const float* d_quadricA,
	const float* d_quadricB,
	const float* d_quadricC,
	const int nbSpheres,
	const float* d_spheresPos,
	const float* d_spheresRadius,
	const float lambda,
	const float tau,
	const float eps,
	const float thetaGap,
	const int distanceMode,
	MembershipEntry* d_outEntries,
	int* d_outCounts)
{

	const int sampleId = blockIdx.x * blockDim.x + threadIdx.x;
	if (sampleId >= nbSamples)
		return;

	const float3 samplePos = load_float3(d_samplesPos, sampleId);
	const float3 sampleProj = load_float3(d_samplesProj, sampleId);
	const float3 sampleNorm = load_float3(d_samplesNorm, sampleId);
	const float* qA = &d_quadricA[16 * sampleId];
	const float* qB = &d_quadricB[4 * sampleId];
	const float qC = d_quadricC[sampleId];
	float bestEnergys[VMAS_MAX_TOPK];
	int bestIndices[VMAS_MAX_TOPK];
	int count = 0;
	for (int s = 0; s < nbSpheres; s++)
	{
		const float3 center = load_float3(d_spheresPos, s); 
		const float radius = d_spheresRadius[s];
		const float energy =
			compute_energy(samplePos, sampleProj, sampleNorm, qA, qB, qC, center, radius, lambda, distanceMode);
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
		return;
	}

	const float Emin = bestEnergys[0];
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

	float denom = 0.f;
	for (int i = 0; i < keepCount; ++i)
	{
		const float w = expf(-bestEnergys[i] / tau);
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
		return;
	}
	if (denom <= 0.f)
		denom = 1.f;
	for (int i = 0; i < keepCount; ++i)
	{
		dst[i].weight /= denom;
	}
	d_outCounts[sampleId] = keepCount;
}

__global__ void null_kernel()
{ /* no-op */
}

extern "C"
cudaError_t cgogn_compute_membership(
	const int nbSamples,
	const float* d_samplesPos,
	const float* d_samplesProj,
	const float* d_samplesNorm,
	const float* d_quadricA,
	const float* d_quadricB,
	const float* d_quadricC,
	const int nbSpheres,
	const float* d_spheresPos,
	const float* d_spheresRadius,
	const float lambda,
	const float tau,
	const float eps,
	const float thetaGap,
	const int distanceMode,
	MembershipEntry* d_outEntries,
	int* d_outCounts)
{
	if (nbSamples == 0 || nbSpheres == 0)
		return cudaSuccess;

	const int blockSize = 256;
	const int gridSize = (nbSamples + blockSize - 1) / blockSize;

	compute_membership_kernel<<<gridSize, blockSize>>>(
		nbSamples,
		d_samplesPos,
		d_samplesProj,
		d_samplesNorm,
		d_quadricA,
		d_quadricB,
		d_quadricC,
		nbSpheres,
		d_spheresPos,
		d_spheresRadius,
		lambda,
		tau,
		eps,
		thetaGap,
		distanceMode,
		d_outEntries,
		d_outCounts);

	return cudaGetLastError();
}