#pragma once

#include <cuda_runtime.h>

#include <cgogn/core/types/cuda_plain_traits.h>

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cgogn
{
namespace cuda
{
template < typename AttributePtr> class CudaAttributePlainBuffer
{
public:
	using Attribute = typename AttributePtr::element_type;
	using ValueType = typename Attribute::value_type;
	using Traits = CudaPlainValueTraits<ValueType>;
	using PlainType = typename Traits::PlainType;

	CudaAttributePlainBuffer() = default;
	explicit CudaAttributePlainBuffer(AttributePtr attribute) : attribute_(std::move(attribute))
	{
		refresh_from_attribute();
	}

	void reset_attributes(AttributePtr attribute, bool refresh = true)
	{
		release_device();
		attribute_ = std::move(attribute);
		host_size_ = 0;
		index_buffer_.clear();

		if (attribute_ && refresh)
			refresh_from_attribute();
	}

	~CudaAttributePlainBuffer()
	{
		release_device();
		release_host();
	}

	CudaAttributePlainBuffer(const CudaAttributePlainBuffer&) = delete;
	CudaAttributePlainBuffer& operator=(const CudaAttributePlainBuffer&) = delete;

	CudaAttributePlainBuffer(CudaAttributePlainBuffer&& other) noexcept
		: attribute_(std::move(other.attribute_)), host_buffer_(other.host_buffer_),
		  host_capacity_(other.host_capacity_), host_size_(other.host_size_),
		  index_buffer_(std::move(other.index_buffer_)), device_ptr_(other.device_ptr_)
	{
		other.host_buffer_ = nullptr;
		other.host_capacity_ = 0;
		other.host_size_ = 0;
		other.device_ptr_ = nullptr;
	}

	CudaAttributePlainBuffer& operator=(CudaAttributePlainBuffer&& other) noexcept
	{

		release_device();
		release_host();

		attribute_ = std::move(other.attribute_);
		host_buffer_ = other.host_buffer_;
		host_capacity_ = other.host_capacity_;
		host_size_ = other.host_size_;
		index_buffer_ = std::move(other.index_buffer_);
		device_ptr_ = other.device_ptr_;

		other.host_buffer_ = nullptr;
		other.host_capacity_ = 0;
		other.host_size_ = 0;
		other.device_ptr_ = nullptr;

		return *this;
	}

	void refresh_from_attribute()
	{
		host_size_ = 0;
		index_buffer_.clear();

		if (!attribute_)
			return;

		const std::size_t required = attribute_->size();
		if (required == 0)
			return;

		reserve_host(required);
		index_buffer_.reserve(required);

		for (auto it = attribute_->begin(), end = attribute_->end(); it != end; ++it)
		{
			index_buffer_.push_back(it.index());
			host_buffer_[host_size_++] = Traits::to_plain(*it);
		}
	}

	void upload(cudaStream_t stream = nullptr)
	{
		release_device();

		if (host_size_ == 0)
			return;

		const std::size_t bytes = host_size_ * sizeof(PlainType);

		cudaError_t st = cudaMalloc(reinterpret_cast<void**>(&device_ptr_), bytes);
		if (st != cudaSuccess)
			throw std::runtime_error(error("cudaMalloc", st));

		st = cudaMemcpyAsync(device_ptr_, host_buffer_, bytes, cudaMemcpyHostToDevice, stream);
		if (st != cudaSuccess)
		{
			cudaFree(device_ptr_);
			device_ptr_ = nullptr;
			throw std::runtime_error(error("cudaMemcpyAsync (H2D)", st));
		}
	}

	void download(cudaStream_t stream = nullptr)
	{
		if (!Traits::can_write_back || device_ptr_ == nullptr || host_size_ == 0)
			return;

		const std::size_t bytes = host_size_ * sizeof(PlainType);

		cudaError_t st = cudaMemcpyAsync(host_buffer_, device_ptr_, bytes, cudaMemcpyDeviceToHost, stream);
		if (st != cudaSuccess)
			throw std::runtime_error(error("cudaMemcpyAsync (D2H)", st));

		// TODO: try to seperate the synchronize call from this function
		//  make sure the copy is finished before writing back to the attribute
		if (stream == nullptr)
		{
			st = cudaDeviceSynchronize();
			if (st != cudaSuccess)
				throw std::runtime_error(error("cudaDeviceSynchronize", st));
		}

		for (std::size_t i = 0; i < index_buffer_.size(); ++i)
			(*attribute_)[index_buffer_[i]] = Traits::from_plain(host_buffer_[i]);
	}

	PlainType* device_data() noexcept
	{
		return device_ptr_;
	}
	const PlainType* device_data() const noexcept
	{
		return device_ptr_;
	}

	std::size_t size() const noexcept
	{
		return host_size_;
	}
	bool empty() const noexcept
	{
		return host_size_ == 0;
	}

	const PlainType* host_data() const noexcept
	{
		return host_buffer_;
	}
	PlainType* host_data() noexcept
	{
		return host_buffer_;
	}

	const std::vector<uint32>& indices() const noexcept
	{
		return index_buffer_;
	}

private:
	void reserve_host(std::size_t n)
	{
		if (n <= host_capacity_)
			return;

		release_host();

		cudaError_t st = cudaMallocHost(reinterpret_cast<void**>(&host_buffer_), n * sizeof(PlainType));
		if (st != cudaSuccess)
			throw std::runtime_error(error("cudaMallocHost", st));

		host_capacity_ = n;
	}

	void release_device()
	{
		if (device_ptr_ != nullptr)
		{
			cudaFree(device_ptr_);
			device_ptr_ = nullptr;
		}
	}

	void release_host()
	{
		if (host_buffer_ != nullptr)
		{
			cudaFreeHost(host_buffer_);
			host_buffer_ = nullptr;
			host_capacity_ = 0;
			host_size_ = 0;
		}
	}

	static std::string error(const char* call, cudaError_t st)
	{
		return std::string(call) + " [" + cudaGetErrorName(st) + "] failed: " + cudaGetErrorString(st);
	}
	AttributePtr attribute_;

	PlainType* host_buffer_ = nullptr; // pinned host memory
	std::size_t host_capacity_ = 0;
	std::size_t host_size_ = 0;

	std::vector<uint32> index_buffer_;
	PlainType* device_ptr_ = nullptr; // device memory
};

} // namespace cuda
} // namespace cgogn
