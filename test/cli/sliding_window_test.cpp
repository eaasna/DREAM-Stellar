#include <fstream>
#include <seqan3/std/ranges>     // range comparisons
#include <string>                // strings
#include <vector>                // vectors

#include "cli_test.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// sliding_window build tests ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEST_P(sliding_window_build, build_with_file)
{
    auto const [number_of_bins, window_size, run_parallel_tmp] = GetParam();
    bool const run_parallel = run_parallel_tmp && number_of_bins >= 32;

    {
        std::ofstream file{"bin_paths.txt"};
	for (size_t i{0}; i < number_of_bins; i++)
        {
	    std::string file_path = cli_test::data("bin_" + std::to_string(i) + ".fasta");
            file << file_path << '\n';
        }
        file << '\n';
    }

    cli_test_result const result = execute_app("sliding_window", "build",
                                                         "--kmer 20",
                                                         "--window ", std::to_string(window_size),
                                                         "--size 100k",
                                                         "--threads ", run_parallel ? "2" : "1",
                                                         "--output index.ibf",
                                                         "bin_paths.txt");
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});

    std::string const expected = string_from_file(ibf_path(number_of_bins, window_size), std::ios::binary);
    std::string const actual = string_from_file("index.ibf", std::ios::binary);

    EXPECT_TRUE(expected == actual);
}


INSTANTIATE_TEST_SUITE_P(build_suite,
                         sliding_window_build,
                         testing::Combine(testing::Values(8), testing::Values(20, 23), testing::Values(true, false)),
                         [] (testing::TestParamInfo<sliding_window_build::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::get<0>(info.param)) + "_bins_" +
                                                std::to_string(std::get<1>(info.param)) + "_window_" +
                                                (std::get<2>(info.param) ? "parallel" : "serial");
                             return name;
                         });

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// sliding_window search tests //////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEST_P(sliding_window_search, search)
{
    auto const [number_of_bins, window_size, number_of_errors, pattern_size, overlap] = GetParam();

    if (window_size == 23 && number_of_errors == 0)
        GTEST_SKIP() << "Needs dynamic threshold correction";

    cli_test_result const result = execute_app("sliding_window", "search",
                                                         "--output search.out",
							 "--pattern", std::to_string(pattern_size),
							 "--overlap", std::to_string(overlap),
                                                         "--error ", std::to_string(number_of_errors),
                                                         "--index ", ibf_path(number_of_bins, window_size),
                                                         "--query ", data("query.fq"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});

    std::string const expected = string_from_file(search_result_path(number_of_bins, window_size, number_of_errors, 
			    pattern_size, overlap), std::ios::binary);
    std::string const actual = string_from_file("search.out");

    EXPECT_EQ(expected, actual);
}

INSTANTIATE_TEST_SUITE_P(search_suite,
                         sliding_window_search,
                         testing::Combine(testing::Values(8), testing::Values(20, 23), testing::Values(0, 1), 
				 testing::Values(50, 100), testing::Values(1, 40)),
                         [] (testing::TestParamInfo<sliding_window_search::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::get<0>(info.param)) + "_bins_" +
                                                std::to_string(std::get<1>(info.param)) + "_window_" +
                                                std::to_string(std::get<2>(info.param)) + "_error_" +
						std::to_string(std::get<3>(info.param)) + "_pattern_" +
						std::to_string(std::get<4>(info.param)) + "_overlap";
                             return name;
                         });
