//  Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <vector>
#include <random>
#include <boost/decimal/decimal128_t.hpp>
#include <boost/decimal/detail/quantize_impl.hpp>
#include <boost/decimal/bid_conversion.hpp>
#include "cuda_managed_ptr.hpp"
#include "stopwatch.hpp"

#include <cuda_runtime.h>

using test_type = boost::decimal::decimal128_t;

__global__ void cuda_test(const test_type *in, const test_type *quantum, test_type *out, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements)
    {
        out[i] = boost::decimal::quantize(in[i], quantum[i]);
    }
}

int main(void)
{
    std::mt19937_64 rng{42};

    int numElements = 50000;
    std::cout << "[Vector operation on " << numElements << " elements]" << std::endl;

    cuda_managed_ptr<test_type> input_vector(numElements);
    cuda_managed_ptr<test_type> quantum_vector(numElements);
    cuda_managed_ptr<test_type> output_vector(numElements);

    // The exponents give results that round, results that gain zeros, and results that do not fit.
    // Half the coefficients have up to 16 digits. A high word below that of 10^34 - 1 gives up to 34 digits.
    std::uniform_int_distribution<std::uint64_t> sig16_dist{0, UINT64_C(9999999999999999)};
    std::uniform_int_distribution<std::uint64_t> high_dist{0, UINT64_C(542101086242751)};
    std::uniform_int_distribution<std::uint64_t> low_dist;
    std::uniform_int_distribution<int> exp_dist{-40, 40};
    for (int i = 0; i < numElements; ++i)
    {
        const boost::int128::uint128_t sig {i % 2 == 0 ? boost::int128::uint128_t {sig16_dist(rng)} : boost::int128::uint128_t {high_dist(rng), low_dist(rng)}};
        input_vector[i] = test_type(sig, exp_dist(rng), i % 3 == 0);
        quantum_vector[i] = test_type(1, exp_dist(rng));
    }

    // Results that do not come from the host quantize: a negative tie, nonzero tails, a coefficient
    // above 64 bits, and a result that gains zeros
    const boost::int128::uint128_t two_words {1U, 1U};
    const boost::int128::uint128_t max_sig {UINT64_C(542101086242752), UINT64_C(4003012203950112767)};
    const boost::int128::uint128_t pow10_32 {UINT64_C(5421010862427), UINT64_C(9632337040368467968)};
    const test_type fixed[][3] {
        {test_type(-125, -2), test_type(1, -1), test_type(-12, -1)},
        {test_type(-135, -2), test_type(1, -1), test_type(-14, -1)},
        {test_type(1251, -3), test_type(1, -1), test_type(13, -1)},
        {test_type(-1249, -3), test_type(1, -1), test_type(-12, -1)},
        {test_type(two_words, -1), test_type(1, 0), test_type(UINT64_C(1844674407370955162), 0)},
        {test_type(max_sig, -2), test_type(1, 0), test_type(pow10_32, 0)},
        {test_type(12, 0), test_type(1, -3), test_type(12000, -3)},
    };
    constexpr int numFixed {static_cast<int>(sizeof(fixed) / sizeof(fixed[0]))};
    for (int i = 0; i < numFixed; ++i)
    {
        input_vector[i] = fixed[i][0];
        quantum_vector[i] = fixed[i][1];
    }

    int threadsPerBlock = 256;
    int blocksPerGrid = (numElements + threadsPerBlock - 1) / threadsPerBlock;
    std::cout << "CUDA kernel launch with " << blocksPerGrid << " blocks of " << threadsPerBlock << " threads" << std::endl;

    watch w;

    cudaDeviceSetLimit(cudaLimitStackSize, 2048);
    cuda_test<<<blocksPerGrid, threadsPerBlock>>>(input_vector.get(), quantum_vector.get(), output_vector.get(), numElements);
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

    std::cout << "CUDA kernel done in: " << w.elapsed() << "s" << std::endl;

    std::vector<test_type> results;
    results.reserve(numElements);
    w.reset();
    for (int i = 0; i < numElements; ++i)
    {
        results.push_back(boost::decimal::quantize(input_vector[i], quantum_vector[i]));
    }
    double t = w.elapsed();

    // The fixed results assume ties to even, the default rounding mode
    const bool ties_to_even {boost::decimal::fegetround() == boost::decimal::rounding_mode::fe_dec_to_nearest};
    for (int i = 0; ties_to_even && i < numFixed; ++i)
    {
        if (boost::decimal::to_bid(output_vector[i]) != boost::decimal::to_bid(fixed[i][2]))
        {
            std::cerr << "Fixed result failed at element " << i << "!" << std::endl;
            return EXIT_FAILURE;
        }
    }

    for (int i = 0; i < numElements; ++i)
    {
        // A NaN result does not compare equal to itself
        const bool same {isnan(results[i]) ? isnan(output_vector[i]) : output_vector[i] == results[i] && samequantumd128(output_vector[i], results[i])};
        if (!same)
        {
            std::cerr << "Result verification failed at element " << i << "!" << std::endl;
            return EXIT_FAILURE;
        }
    }

    std::cout << "Test PASSED, normal calculation time: " << t << "s" << std::endl;
    std::cout << "Done\n";

    return 0;
}
