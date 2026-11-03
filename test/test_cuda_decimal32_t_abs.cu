//  Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <vector>
#include <random>
#include <limits>
#include <boost/decimal/decimal32_t.hpp>
#include <boost/decimal/bid_conversion.hpp>
#include "cuda_managed_ptr.hpp"

#include <cuda_runtime.h>

using test_type = boost::decimal::decimal32_t;

__global__ void cuda_test(const test_type *in, test_type *out, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements)
    {
        out[i] = boost::decimal::abs(in[i]);
    }
}

int main(void)
{
    std::mt19937_64 rng{42};

    int numElements = 50000;
    std::cout << "[Vector operation on " << numElements << " elements]" << std::endl;

    cuda_managed_ptr<test_type> input_vector(numElements);
    cuda_managed_ptr<test_type> output_vector(numElements);

    // Most inputs are the negation of the expected result, so the kernel clears a sign at run time.
    // Each fourth input is the expected result, so a positive value and a quiet NaN do not change.
    std::vector<test_type> expected(numElements);
    std::uniform_int_distribution<std::int32_t> sig_dist{0, 9999999};
    std::uniform_int_distribution<int> exp_dist{-90, 90};
    for (int i = 0; i < numElements; ++i)
    {
        expected[i] = test_type(sig_dist(rng), exp_dist(rng));
    }
    expected[0] = test_type(0);
    expected[1] = std::numeric_limits<test_type>::infinity();
    expected[2] = std::numeric_limits<test_type>::denorm_min();
    expected[3] = std::numeric_limits<test_type>::quiet_NaN();
    for (int i = 0; i < numElements; ++i)
    {
        input_vector[i] = i % 4 == 3 ? expected[i] : -expected[i];
    }

    int threadsPerBlock = 256;
    int blocksPerGrid = (numElements + threadsPerBlock - 1) / threadsPerBlock;
    std::cout << "CUDA kernel launch with " << blocksPerGrid << " blocks of " << threadsPerBlock << " threads" << std::endl;

    cuda_test<<<blocksPerGrid, threadsPerBlock>>>(input_vector.get(), output_vector.get(), numElements);
    cudaError_t err {cudaGetLastError()};
    if (err == cudaSuccess)
    {
        err = cudaDeviceSynchronize();
    }
    if (err != cudaSuccess)
    {
        std::cerr << "CUDA error: " << cudaGetErrorString(err) << std::endl;
        return EXIT_FAILURE;
    }

    for (int i = 0; i < numElements; ++i)
    {
        if (boost::decimal::to_bid(output_vector[i]) != boost::decimal::to_bid(expected[i]))
        {
            std::cerr << "Result verification failed at element " << i << "!" << std::endl;
            return EXIT_FAILURE;
        }
    }

    std::cout << "Test PASSED" << std::endl;
    std::cout << "Done\n";

    return 0;
}
