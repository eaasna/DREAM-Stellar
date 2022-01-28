#include <valik/search/query_record.hpp>
#include <valik/search/search_setup.hpp>

namespace valik::app
{

//-----------------------------
//
// Setup IBF and launch multithreaded search.
//
//-----------------------------
template <seqan3::data_layout ibf_data_layout>
void run_program(search_arguments const &arguments, search_time_statistics & time_statistics)
{
    seqan3::interleaved_bloom_filter<ibf_data_layout> ibf{};

    auto cereal_worker = [&]()
    {
        auto start = std::chrono::high_resolution_clock::now();
        load_ibf(ibf, arguments.ibf_file);
        auto end = std::chrono::high_resolution_clock::now();
        time_statistics.ibf_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    };

    auto cereal_handle = std::async(std::launch::async, cereal_worker);

    seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin{arguments.query_file};
    std::vector<query_record> query_records{};

    threshold const threshold_data = make_threshold_data(arguments);

    sync_out synced_out{arguments.out_file};

    size_t record_id = 0u;
    for (auto &&chunked_records : fin | seqan3::views::chunk((1ULL << 20) * 10))
    {
        query_records.clear();
        auto start = std::chrono::high_resolution_clock::now();
        for (auto && query_record: chunked_records)
        {
            query_records.emplace_back(record_id, std::move(query_record.id()), std::move(query_record.sequence()));
            ++record_id;
        }
        auto end = std::chrono::high_resolution_clock::now();
        time_statistics.reads_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

        cereal_handle.wait();

        start = std::chrono::high_resolution_clock::now();
        write_output_file_parallel(ibf, arguments, query_records, threshold_data, synced_out);
        end = std::chrono::high_resolution_clock::now();
        time_statistics.compute_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    }
}

void valik_search(search_arguments const & arguments)
{
    search_time_statistics time_statistics{};

    if (arguments.compressed)
        run_program<seqan3::data_layout::compressed>(arguments, time_statistics);
    else
        run_program<seqan3::data_layout::uncompressed>(arguments, time_statistics);

    if (arguments.write_time)
        write_time_statistics(time_statistics, arguments);

    return;
}

} // namespace valik::app
