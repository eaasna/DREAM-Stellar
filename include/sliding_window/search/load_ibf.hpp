#pragma once

#include <chrono>
#include <seqan3/std/filesystem>

#include <sliding_window/shared.hpp>

namespace sliding_window
{

template <typename t>
void load_ibf(t & ibf, search_arguments const & arguments, double & ibf_io_time)
{
    static uint8_t kmer_size{};
    static uint32_t window_size{};

    std::ifstream is{arguments.ibf_file, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};

    auto start = std::chrono::high_resolution_clock::now();
    iarchive(kmer_size);
    iarchive(window_size);
    iarchive(ibf);
    auto end = std::chrono::high_resolution_clock::now();

    ibf_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
}

} // namespace sliding_window
