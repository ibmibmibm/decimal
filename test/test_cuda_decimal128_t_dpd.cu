//  Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <vector>
#include <random>
#include <boost/decimal/decimal128_t.hpp>
#include <boost/decimal/dpd_conversion.hpp>
#include "cuda_managed_ptr.hpp"
#include "stopwatch.hpp"

#include <cuda_runtime.h>

using test_type = boost::decimal::decimal128_t;

using bits_type = boost::int128::uint128_t;

// The device computes each declet, so the host tables check its bits
__global__ void cuda_test(const test_type *in, test_type *out, bits_type *out_bits, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements)
    {
        auto bits = boost::decimal::to_dpd(in[i]);
        out_bits[i] = bits;
        out[i] = boost::decimal::from_dpd<test_type>(bits);
    }
}

int main(void)
{
    std::mt19937_64 rng{42};

    int numElements = 50000;
    std::cout << "[Vector operation on " << numElements << " elements]" << std::endl;

    cuda_managed_ptr<test_type> input_vector(numElements);
    cuda_managed_ptr<test_type> output_vector(numElements);
    cuda_managed_ptr<bits_type> output_bits(numElements);

    // Coefficients up to 34 digits, and a signaling NaN with a payload
    const bits_type limit {UINT64_C(542101086242752), UINT64_C(4003012203950112768)};
    std::uniform_int_distribution<int> exp_dist{-6176, 6111};
    for (int i = 0; i < numElements; ++i)
    {
        input_vector[i] = test_type(bits_type{rng(), rng()} % limit, exp_dist(rng));
    }
    input_vector[0] = boost::decimal::from_bid<test_type>(bits_type{UINT64_C(0xFE00000000000000), UINT64_C(123)});

    int threadsPerBlock = 256;
    int blocksPerGrid = (numElements + threadsPerBlock - 1) / threadsPerBlock;
    std::cout << "CUDA kernel launch with " << blocksPerGrid << " blocks of " << threadsPerBlock << " threads" << std::endl;

    watch w;

    cudaDeviceSetLimit(cudaLimitStackSize, 2048);
    cuda_test<<<blocksPerGrid, threadsPerBlock>>>(input_vector.get(), output_vector.get(), output_bits.get(), numElements);
    cudaDeviceSynchronize();

    std::cout << "CUDA kernel done in: " << w.elapsed() << "s" << std::endl;

    std::vector<test_type> results;
    std::vector<bits_type> results_bits;
    results.reserve(numElements);
    results_bits.reserve(numElements);
    w.reset();
    for (int i = 0; i < numElements; ++i)
    {
        auto bits = boost::decimal::to_dpd(input_vector[i]);
        results_bits.push_back(bits);
        results.push_back(boost::decimal::from_dpd<test_type>(bits));
    }
    double t = w.elapsed();

    for (int i = 0; i < numElements; ++i)
    {
        if (output_bits[i] != results_bits[i] ||
            boost::decimal::to_bid(output_vector[i]) != boost::decimal::to_bid(results[i]))
        {
            std::cerr << "Result verification failed at element " << i << "!" << std::endl;
            return EXIT_FAILURE;
        }
    }

    std::cout << "Test PASSED, normal calculation time: " << t << "s" << std::endl;
    std::cout << "Done\n";

    return 0;
}
