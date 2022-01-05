#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/slice.hpp> // provides views::slice

#include <sliding_window/search/compute_simple_model.hpp>
#include <sliding_window/search/write_output_file_parallel.hpp>
#include <sliding_window/search/load_ibf.hpp>
#include <sliding_window/search/query_result.hpp>
#include <sliding_window/search/sync_out.hpp>

#include <indexed_minimiser_hash.hpp>

namespace sliding_window
{

//-----------------------------
//
// Precalculate vector with all pattern begin positions in read.
// This adds a small computational overhead but makes the code much more readable.
//
//-----------------------------
std::vector<size_t> precalculate_begin(size_t const read_len, uint64_t const pattern_size, uint64_t const overlap)
{
    std::vector<size_t> begin_vector;
    for (size_t i = 0; i <= read_len - pattern_size;
            i = i + pattern_size - overlap)
    {
        begin_vector.push_back(i);
    }
    if (begin_vector.back() < read_len - pattern_size)
    {
        // last pattern might have a smaller overlap to make sure the end of the read is covered
        begin_vector.push_back(read_len - pattern_size);
    }

    return begin_vector;
}

// TODO: make these generic for both build and search
struct search_time_statistics
{
    double ibf_io_time{0.0};
    double reads_io_time{0.0};
    double compute_time{0.0};
};

void write_time_statistics(search_time_statistics const & time_statistics, search_arguments const & arguments)
{
    std::filesystem::path file_path{arguments.out_file};
    file_path += ".time";
    std::ofstream file_handle{file_path};
    file_handle << "IBF I/O\tReads I/O\tCompute\n";
    file_handle << std::fixed
                << std::setprecision(2)
                << time_statistics.ibf_io_time << '\t'
                << time_statistics.reads_io_time << '\t'
                << time_statistics.compute_time;

}

//-----------------------------
//
// Position of a sliding window (i.e pattern) on the read and threshold for local match.
//
//-----------------------------
struct pattern_bounds
{
    size_t first_index;
    size_t last_index;
    size_t threshold;
};

template <typename span_vec_t>
pattern_bounds make_pattern_bounds(size_t const & begin, search_arguments const & arguments, span_vec_t const & window_span_begin, 
                                                            span_vec_t const & window_span_end, threshold const & threshold_data)
{
    auto pattern = pattern_bounds{};

    // indices for the first and last minimiser of the current sliding window
    // std::vector<size_t>::iterator lower_it, upper_it;
    auto lower_it = std::lower_bound(window_span_begin.begin(), window_span_begin.end(), begin);
    auto upper_it = std::upper_bound(window_span_end.begin(), window_span_end.end(),
                                begin + arguments.pattern_size - 1); // - 1 because of 0 based indexing

    pattern.first_index = lower_it - window_span_begin.begin();
    pattern.last_index = upper_it - window_span_end.begin() - 1; // - 1 because the upper bound returns the first el that is greater

    size_t const minimiser_count = pattern.last_index - pattern.first_index + 1;

    pattern.threshold = arguments.treshold_was_set ? 
                            static_cast<size_t>(minimiser_count * arguments.threshold) : threshold_data.kmers_per_window == 1 ? 
                            threshold_data.kmer_lemma : threshold_data.precomp_thresholds[std::min(minimiser_count < threshold_data.min_number_of_minimisers ? 
                            0 : minimiser_count - threshold_data.min_number_of_minimisers,  
                            threshold_data.max_number_of_minimisers - threshold_data.min_number_of_minimisers)] + 2; // enrico recommended decreasing this value

    return pattern;
}

//-----------------------------
//
// Count matching k-mers and return matching bins for one pattern.
//
//-----------------------------
template <typename binning_bitvector_t>
std::set<size_t> find_pattern_bins(pattern_bounds const & pattern, size_t const & bin_count, binning_bitvector_t const & counting_table)
{
    std::set<size_t> pattern_hits{};    // bin hits for one pattern i.e sliding window

    // counting vector for the current pattern
    seqan3::counting_vector<uint8_t> total_counts(bin_count, 0);

    for (size_t i = pattern.first_index; i <= pattern.last_index; i++)
        total_counts += counting_table[i];
    for (size_t current_bin = 0; current_bin < total_counts.size(); current_bin++)
    {
        auto &&count = total_counts[current_bin];
        if (count >= pattern.threshold)
        {
            // the result_set is a union of results from all sliding windows of a read
            pattern_hits.insert(current_bin);
        }
    }

    return pattern_hits;
}

//-----------------------------
//
// Search a batch of reads (start; end) in the IBF
//
// TODO: pass slice of records directly instead of start, end and records
//-----------------------------
template <typename ibf_t, typename rec_vec_t>
auto worker(size_t const start, size_t const end, ibf_t const & ibf, search_arguments const & arguments, 
                                                                    rec_vec_t const & records, threshold const & threshold_data)
{
    // concurrent invocations of the membership agent are not thread safe
    // agent has to be created for each thread
    auto &&agent = ibf.membership_agent();
    size_t bin_count = ibf.bin_count();

    std::vector<query_result> thread_result{}; // set of query results processed by one thread
    std::set<size_t> sequence_hits{};          // bin hits for one sequence

    // vector holding all the minimisers and their starting position for the read
    using minimiser_vec_t = std::vector<std::tuple<uint64_t, size_t>>;
    minimiser_vec_t minimiser;

    auto hash_tuple_view = indexed_minimiser_hash(seqan3::ungapped{arguments.kmer_size},
                                                    window_size{arguments.window_size},
                                                    seed{adjust_seed(arguments.kmer_size)});

    for (auto &&[id, seq] : records | seqan3::views::slice(start, end))
    {
        sequence_hits.clear();

        minimiser = seq | hash_tuple_view | seqan3::views::to<minimiser_vec_t>;

        //-----------------------------
        //
        // For each read the begin_vector shows the beginning of each sliding window (pattern)
        //
        // If 	read_len = 150
        // 	pattern_size = 50
        // 	overlap = 30
        //
        // 	begin_vector = {0, 30, 60, 90, 100}
        //
        //-----------------------------
        auto begin_vector = precalculate_begin(seq.size(), arguments.pattern_size, arguments.overlap);

        //-----------------------------
        //
        // Table of counting vectors newly created for each read
        //	rows: each minimiser of read
        // 	columns: each bin of IBF
        //
        //-----------------------------
        using binning_bitvector_t = typename std::remove_cvref_t<decltype(ibf)>::membership_agent::binning_bitvector;
        std::vector<binning_bitvector_t> counting_table(minimiser.size(),
                                                        binning_bitvector_t(bin_count));

        // the beginning of the first window this minimiser is in
        std::vector<size_t> window_span_begin(minimiser.size(), 0);
        // the end of the last window this minimiser is in
        std::vector<size_t> window_span_end(minimiser.size(), 0);

        for (size_t i{1}; i < minimiser.size(); i++)
        {
            const auto &[min, start_pos] = minimiser[i];
            window_span_begin[i] = start_pos;
            size_t end_pos = start_pos + arguments.window_size - 2;
            window_span_end[i - 1] = end_pos;
            counting_table[i].raw_data() |= agent.bulk_contains(min).raw_data();
        }

        window_span_end[minimiser.size() - 1] = seq.size() - 1;
        minimiser.clear();

        //-----------------------------
        //
        // If seq = CGCAAAACGCGGC
        // 	p = 12
        // 	w = 8
        // 	k = 4
        //
        // minimiser 		= (AAAA; 3), (AAAC; 4), (AACG, 5)
        // window_span_begin 	=     0		4	   5
        // window_span_end 	=     10	11	   12
        //
        // minimiser	span
        // AAAA		CGCAAAACGCG
        // AAAC		AAACGCGG
        // AACG		AACGCGGC
        //
        //-----------------------------
        for (auto begin : begin_vector)
        {
            auto const pattern = make_pattern_bounds(begin, arguments, window_span_begin, window_span_end, threshold_data);
            auto const pattern_hits = find_pattern_bins(pattern, bin_count, counting_table);
            sequence_hits.insert(pattern_hits.begin(), pattern_hits.end());
        }

        query_result query_result(id, sequence_hits);
        thread_result.emplace_back(query_result);
    }

    return thread_result;
};

} // namespace sliding_window
