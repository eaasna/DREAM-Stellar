#pragma once

#include <chrono>
#include <future>
#include <vector>

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <valik/shared.hpp>
#include <valik/search/local_prefilter.hpp>
#include <valik/search/query_record.hpp>
#include <valik/search/query_result.hpp>
#include <valik/search/sync_out.hpp>

#include <raptor/threshold/threshold.hpp>

namespace valik::app
{

template <seqan3::data_layout ibf_data_layout>
inline void write_output_file_parallel(seqan3::interleaved_bloom_filter<ibf_data_layout> const & ibf,
                                       search_arguments const & arguments,
                                       std::vector<query_record> const & records,
                                       raptor::threshold::threshold const & thresholder,
                                       sync_out & synced_out)
{
    using query_list = std::vector<std::pair<std::string, std::vector<seqan3::dna4>>>;
    using task_future_t = std::future<std::pair<std::vector<query_result>, std::map<size_t, query_list>>>;
    static_assert(std::same_as<task_future_t,
                               decltype(std::async(std::launch::async, local_prefilter, std::span<query_record const>{}, ibf, arguments, thresholder))>);

    std::vector<task_future_t> tasks;
    size_t const num_records = records.size();
    size_t const records_per_thread = num_records / arguments.threads;

    for (size_t i = 0; i < arguments.threads; ++i)
    {
        size_t const start = records_per_thread * i;
        size_t const end = std::min(start + records_per_thread, num_records);
        std::span<query_record const> records_slice{&records[start], &records[end]};

        // The following calls `local_prefilter(records, ibf, arguments, threshold)` on a thread.
        // Note: local_prefilter is a function object which is created from the local_function_fn struct
        tasks.emplace_back(std::async(std::launch::async, local_prefilter, records_slice, ibf, arguments, thresholder));
    }


    std::filesystem::create_directories("matches");
    for (task_future_t & task : tasks)
    {
        std::string result_string{};
        std::pair<std::vector<query_result>, std::map<size_t, query_list>> const & thread_result = task.get();
        auto const & thread_query_result = thread_result.first;
        auto const & thread_bin_result = thread_result.second;

        for (query_result const & query_result : thread_query_result)
        {
            result_string.clear();
            result_string += query_result.get_id();
            result_string += '\t';

            for (size_t const bin : query_result.get_hits())
            {
                result_string += std::to_string(bin);
                result_string += ',';
            }

            result_string += '\n';
            synced_out.write(result_string);
        }

        // TODO: open multiple bin-query synced_out to write all bin-query to file
        for (auto & [bin_id, matches] : thread_bin_result)
        {
            std::string result_string{};
            std::filesystem::path bin_query_path = "matches/" + std::to_string(bin_id);
            sync_out synced_out{bin_query_path};

            for (auto & [id, seq] : matches)
            {
                result_string += ">";
                result_string += id;
                result_string += '\n';
                for (auto & c : seq)
                {
                    result_string += c.to_char();
                }
                result_string += '\n';
                synced_out.write(result_string);
                result_string.clear();
            }
        }

    }
}

} // namespace valik::app
