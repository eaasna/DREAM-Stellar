#include <filesystem>
#include <fstream>
#include <ranges>     // range comparisons
#include <string>                // strings
#include <vector>                // vectors

#include "cli_test.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// DREAM short shared memory ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEST_P(dream_short_search, short_shared_mem)
{
    auto const [number_of_bins, window_size, number_of_errors] = GetParam();
    size_t pattern_size = 50;
    float error_rate = (float) number_of_errors / (float) pattern_size;

    setup_tmp_dir();
    setenv("VALIK_MERGE", "cat", true);

    std::filesystem::path ref_meta_path = data("ref_meta.txt");
    valik::database_metadata reference(ref_meta_path, false);
    std::filesystem::path seg_meta_path = data("seg_meta150overlap" + std::to_string(number_of_bins) + "bins.txt");
    std::filesystem::path index_path = ibf_path(number_of_bins, window_size);

    cli_test_result const build = execute_app("valik", "build",
                                                        data("ref.fasta"),
                                                        "--window", std::to_string(window_size),
                                                        "--kmer 13",
                                                        "--size 32k",
                                                        "--from-segments",
                                                        "--ref-meta", ref_meta_path,
                                                        "--seg-meta", seg_meta_path,
                                                        "--output ", index_path);

    cli_test_result const result = execute_app("valik", "search",
                                                        "--output search.gff",
                                                        "--pattern ", std::to_string(pattern_size),
                                                        "--overlap ", std::to_string(pattern_size - 1),
                                                        "--error-rate ", std::to_string(error_rate),
                                                        "--index ", index_path,
                                                        "--query ", data("query.fastq"),
                                                        "--ref-meta", ref_meta_path,
                                                        "--seg-meta", seg_meta_path);
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{"Launching stellar search on a shared memory machine...\nLoaded 3 database sequences.\n"});
    EXPECT_EQ(result.err, std::string{"WARNING: Non-unique query ids. Output can be ambiguous.\n"});

    auto distributed = valik::read_stellar_output(search_result_path(number_of_bins, window_size, number_of_errors), reference, std::ios::binary);
    auto local = valik::read_stellar_output("search.gff", reference);

    compare_gff_out(distributed, local);
}

INSTANTIATE_TEST_SUITE_P(short_shared_memory_dream_suite,
                         dream_short_search,
                         testing::Combine(testing::Values(4, 16), testing::Values(13, 15), testing::Values(1)),
                         [] (testing::TestParamInfo<dream_short_search::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::get<0>(info.param)) + "_bins_" +
                                                std::to_string(std::get<1>(info.param)) + "_window_" +
                                                std::to_string(std::get<2>(info.param)) + "_error";
                             return name;
                         });

TEST_F(dream_short_search, no_matches)
{
    setup_tmp_dir();
    setenv("VALIK_MERGE", "cat", true);

    size_t number_of_bins = 4;
    size_t window_size = 15;
    cli_test_result const build = execute_app("valik", "build",
                                                        data("ref.fasta"),
                                                        "--window ", std::to_string(window_size),
                                                        "--kmer 13",
                                                        "--size 32k",
                                                        "--from-segments",
                                                        "--ref-meta", data("ref_meta.txt"),
                                                        "--seg-meta", data("seg_meta150overlap" + std::to_string(number_of_bins) + "bins.txt"),
                                                        "--output ", ibf_path(number_of_bins, window_size));

    cli_test_result const result = execute_app("valik", "search",
                                                        "--output search.gff",
                                                        "--pattern 50",
                                                        "--overlap 49",
                                                        "--error-rate 0",
                                                        "--index ", ibf_path(number_of_bins, window_size),
                                                        "--query ", data("dummy_reads.fastq"),
                                                        "--ref-meta", data("ref_meta.txt"),
                                                        "--seg-meta", data("seg_meta150overlap" + std::to_string(number_of_bins) + "bins.txt"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{"Launching stellar search on a shared memory machine...\nLoaded 3 database sequences.\n"});
    EXPECT_EQ(result.err, std::string{}); // Stellar shortens query IDs until the first whitespace

    auto actual = string_list_from_file("search.gff");

    EXPECT_EQ(0, actual.size());
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// DREAM split shared memory ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEST_P(dream_split_search, split_shared_mem)
{
    auto const [number_of_bins, window_size, number_of_errors] = GetParam();

    setup_tmp_dir();
    setenv("VALIK_MERGE", "cat", true);
    size_t pattern_size = 50;
    float error_rate = (float) number_of_errors / (float) pattern_size;

    size_t query_seg_count = 60;
    std::filesystem::path ref_meta_path = data("ref_meta.txt");
    valik::database_metadata reference(ref_meta_path, false);
    std::filesystem::path seg_meta_path = data("seg_meta150overlap" + std::to_string(number_of_bins) + "bins.txt");
    std::filesystem::path index_path = ibf_path(number_of_bins, window_size);
    std::filesystem::path query_meta_path = data("query_meta.txt");
    std::filesystem::path query_seg_meta_path = data("query_seg_meta.txt");

    cli_test_result const split_query = execute_app("valik", "split",
                                                        data("query.fastq"),
                                                        "--db-meta", query_meta_path,
                                                        "--seg-meta", query_seg_meta_path,
                                                        "--seg-count ", std::to_string(query_seg_count),
                                                        "--overlap 0");


    cli_test_result const build = execute_app("valik", "build",
                                                        data("ref.fasta"),
                                                        "--window", std::to_string(window_size),
                                                        "--kmer 13",
                                                        "--size 32k",
                                                        "--from-segments",
                                                        "--ref-meta", ref_meta_path,
                                                        "--seg-meta", seg_meta_path,
                                                        "--output ", index_path);

    cli_test_result const search = execute_app("valik", "search",
                                                        "--output search.gff",
                                                        "--pattern ", std::to_string(pattern_size),
                                                        "--overlap ", std::to_string(pattern_size - 1),
                                                        "--error-rate ", std::to_string(error_rate),
                                                        "--index ", index_path,
                                                        "--query ", data("query.fastq"),
                                                        "--ref-meta", ref_meta_path,
                                                        "--seg-meta", seg_meta_path,
                                                        "--query-meta", query_seg_meta_path);

    EXPECT_EQ(search.exit_code, 0);
    EXPECT_EQ(search.out, std::string{"Launching stellar search on a shared memory machine...\nLoaded 3 database sequences.\n"});
    EXPECT_EQ(search.err, std::string{"WARNING: Non-unique query ids. Output can be ambiguous.\n"});

    auto distributed = valik::read_stellar_output(search_result_path(number_of_bins, window_size, number_of_errors), reference, std::ios::binary);
    auto local = valik::read_stellar_output("search.gff", reference);

    compare_gff_out(distributed, local);
}

INSTANTIATE_TEST_SUITE_P(split_shared_memory_dream_suite,
                         dream_split_search,
                         testing::Combine(testing::Values(4, 16), testing::Values(13, 15), testing::Values(1)),
                         [] (testing::TestParamInfo<dream_split_search::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::get<0>(info.param)) + "_bins_" +
                                                std::to_string(std::get<1>(info.param)) + "_window_" +
                                                std::to_string(std::get<2>(info.param)) + "_error";
                             return name;
                         });
