// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <iostream>
#include <cstring>
#include <boost/decimal/decimal64_t.hpp>
#include <boost/decimal/charconv.hpp>
#include "cuda_managed_ptr.hpp"

#include <cuda_runtime.h>

using test_type = boost::decimal::decimal64_t;

// Each result has a slot with one guard byte before the output range and guard bytes after it
constexpr int slot_size {48};

__global__ void cuda_test(const test_type *in, const int *sizes, char *out, int *lengths, int *errors, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements)
    {
        char* first {out + i * slot_size + 1};
        const auto r {boost::decimal::to_chars(first, first + sizes[i], in[i])};
        lengths[i] = static_cast<int>(r.ptr - first);
        errors[i] = static_cast<int>(r.ec);
    }
}

struct test_case
{
    test_type value;
    int size;
    const char* text; // nullptr for the error value_too_large
};

int main(void)
{
    // The exact text of the general format, as the round trip of the other tests does not show the notation
    const test_case cases[] {
        // One group of eight digits with leading zeros, two groups, a sign, and a value above 32 bits
        {test_type(100000001), 39, "100000001"},
        {test_type(UINT64_C(1234567890123456)), 39, "1234567890123456"},
        {test_type(-100000001), 39, "-100000001"},
        {test_type(UINT64_C(4294967296)), 39, "4294967296"},
        {test_type(123456789, -4), 39, "12345.6789"},
        // Ranges of the exact length
        {test_type(100000001), 9, "100000001"},
        {test_type(-100000001), 10, "-100000001"},
        // Ranges one character too short
        {test_type(100000001), 8, nullptr},
        {test_type(-100000001), 9, nullptr},
    };
    constexpr int numElements {static_cast<int>(sizeof(cases) / sizeof(cases[0]))};

    cuda_managed_ptr<test_type> input_vector(numElements);
    cuda_managed_ptr<int> size_vector(numElements);
    cuda_managed_ptr<char> output_vector(numElements * slot_size);
    cuda_managed_ptr<int> length_vector(numElements);
    cuda_managed_ptr<int> error_vector(numElements);

    for (int i = 0; i < numElements; ++i)
    {
        input_vector[i] = cases[i].value;
        size_vector[i] = cases[i].size;
    }
    std::memset(output_vector.get(), '#', static_cast<std::size_t>(numElements * slot_size));

    cuda_test<<<1, numElements>>>(input_vector.get(), size_vector.get(), output_vector.get(), length_vector.get(), error_vector.get(), numElements);
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
        const char* slot {output_vector.get() + i * slot_size};
        bool guards_kept {slot[0] == '#'};
        for (int j = 1 + cases[i].size; j < slot_size; ++j)
        {
            guards_kept = guards_kept && slot[j] == '#';
        }

        bool same {};
        if (cases[i].text == nullptr)
        {
            same = error_vector[i] == static_cast<int>(std::errc::value_too_large) && length_vector[i] == cases[i].size;
        }
        else
        {
            const auto length {static_cast<int>(std::strlen(cases[i].text))};
            same = error_vector[i] == 0 && length_vector[i] == length && std::memcmp(slot + 1, cases[i].text, static_cast<std::size_t>(length)) == 0;
        }

        if (!same || !guards_kept)
        {
            std::cerr << "Result verification failed at element " << i << "!" << std::endl;
            return EXIT_FAILURE;
        }
    }

    std::cout << "Test PASSED" << std::endl;
    std::cout << "Done\n";

    return 0;
}
