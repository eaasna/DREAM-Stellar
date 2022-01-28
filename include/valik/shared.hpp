#pragma once

#include <seqan3/std/filesystem>
#include <vector>

#include <seqan3/io/sequence_file/input.hpp>

namespace valik
{

inline constexpr static uint64_t adjust_seed(uint8_t const kmer_size, uint64_t const seed = 0x8F3F73B5CF1C9ADEULL) noexcept
{
    return seed >> (64u - 2u * kmer_size);
}

//!\brief Strong type for passing the window size.
struct window { uint64_t v; };
//!\brief Strong type for passing the kmer size.
struct kmer { uint8_t v; };
//!\brief Strong type for passing number of bins.
struct bins { uint64_t v; };
//!\brief Strong type for passing number of bits.
struct bits { uint64_t v; };
//!\brief Strong type for passing number of hash functions.
struct hashes { uint64_t v; };

struct dna4_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

struct split_arguments
{
    std::filesystem::path ref_file{};
    std::filesystem::path ref_out{"reference_metadata.txt"};
    std::filesystem::path seg_out{"reference_segments.txt"};

    size_t overlap{150};
    size_t min_bins{64};
    size_t min_len{1000000};
};

struct build_arguments
{
    uint32_t window_size{23u};
    uint8_t kmer_size{20u};
    uint8_t threads{1u};

    std::vector<std::vector<std::filesystem::path>> bin_path{};
    std::filesystem::path bin_file{};
    std::filesystem::path out_path{"./"};
    std::string size{};
    uint64_t bins{64};
    uint64_t bits{4096};
    uint64_t hash{2};
    bool compressed{false};
    
    bool from_segments{false};
    std::filesystem::path seg_path{};
    std::filesystem::path ref_meta_path{};
};

struct search_arguments
{
    uint32_t window_size{23u};
    uint8_t kmer_size{20u};
    uint8_t threads{1u};

    std::filesystem::path query_file{};
    std::filesystem::path ibf_file{};
    std::filesystem::path out_file{"search.out"};
    double tau{0.99};
    double threshold{};
    bool treshold_was_set{false};
    uint64_t pattern_size{};
    uint64_t overlap{};
    uint8_t errors{0};
    bool compressed{false};
    bool write_time{false};
};

} // namespace valik
