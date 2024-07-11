#pragma once

#include <utilities/alphabet_wrapper/matcher/stellar_matcher.hpp>
#include <utilities/alphabet_wrapper/seqan/alphabet.hpp>
#include <valik/search/producer_threads_parallel.hpp>
#include <valik/search/search_time_statistics.hpp>
#include <valik/shared.hpp>

#include <stellar/utils/stellar_app_runtime.hpp>
#include <stellar/io/import_sequence.hpp>

#include <seqan/seq_io.h>

#include <seqan3/io/sequence_file/all.hpp>

namespace valik::app
{

/**
 * @brief Function that sends chunks of queries to the prefilter which then writes shopping carts onto disk.
 *
 * @param arguments Command line arguments.
 * @param ibf Interleaved Bloom Filter.
 * @param thresholder Threshold for number of shared k-mers.
 * @param queue Shopping cart queue for load balancing between prefiltering and Stellar search.
 */
template <typename ibf_t, typename cart_queue_t>
void iterate_distributed_queries(search_arguments const & arguments,
                                 ibf_t const & ibf,
                                 raptor::threshold::threshold const & thresholder,
                                 cart_queue_t & queue)
{
    using fields = seqan3::fields<seqan3::field::id, seqan3::field::seq>;
    std::vector<query_record> query_records{};
    seqan3::sequence_file_input<dna4_traits, fields> fin{arguments.query_file};
    for (auto &&chunked_records : fin | seqan3::views::chunk((1ULL << 20) * 10))
    {
        query_records.clear();
        for (auto && fasta_record: chunked_records)
            query_records.emplace_back(std::move(fasta_record.id()), std::move(fasta_record.sequence()));

        prefilter_queries_parallel(ibf, arguments, query_records, thresholder, queue);
    }
}

/**
 * @brief Function that creates a query record from each query sequence and sends it to Stellar search.
 *
 * @param ref_seg_count Number of reference segments i.e the distribution granularity.
 * @param arguments Command line arguments.
 * @param queue Shopping cart queue for sending queries over to Stellar search.
 */
template <typename TSequence>
void iterate_all_queries(size_t const ref_seg_count,
                         search_arguments const & arguments,
                         cart_queue<shared_query_record<TSequence>> & queue)
{
    using TId = seqan2::CharString;
    std::vector<shared_query_record<TSequence>> query_records{};
    
    constexpr uint64_t chunk_size = (1ULL << 20) * 10;

    seqan2::SeqFileIn inSeqs;
    if (!open(inSeqs, arguments.query_file.c_str()))
    {
        throw std::runtime_error("Failed to open " + arguments.query_file.string() + " file.");
    }

    std::set<TId> uniqueIds; // set of short IDs (cut at first whitespace)
    bool idsUnique = true;

    size_t seqCount{0};

    {   // try stellar_matcher
        seqan3::sequence_file_input<dna4_adaptor_traits> fin{std::istringstream{arguments.query_file}, seqan3::format_fasta{}};
        using record_type = typename decltype(fin)::record_type;
        using sequence_type = std::remove_cvref_t<decltype(std::declval<record_type>().sequence())>;

        std::vector<record_type> rec_vec{};
        for (auto & record : fin)
        {
            rec_vec.emplace_back(record);
        }

        std::vector<sequence_type> seq_vec(rec_vec.size());
        for (auto seq : rec_vec | std::views::transform([](record_type record) { return record.sequence(); }))
        {
            seq_vec.emplace_back(seq);
        }
        jst::contrib::stellar_matcher<sequence_type> matcher(seq_vec, (double) arguments.error_rate, (unsigned) arguments.minLength);
    }
    
    for (; !atEnd(inSeqs); ++seqCount)
    {
        TSequence seq{};
        TId id{};
        readRecord(id, seq, inSeqs);
        idsUnique &= stellar::_checkUniqueId(uniqueIds, id);

        query_records.emplace_back(std::move(seq), seqan2::toCString(std::move(id)));

        if (query_records.size() > chunk_size)
        {
            search_all_parallel<shared_query_record<TSequence>>(ref_seg_count, arguments, query_records, queue);
            query_records.clear();
        }
    }

    if (!idsUnique)
        std::cerr << "WARNING: Non-unique query ids. Output can be ambiguous.\n";

    search_all_parallel<shared_query_record<TSequence>>(ref_seg_count, arguments, query_records, queue);    
}

/**
 * @brief Function that creates short query records from fasta file input and sends them for prefiltering.
 *
 * @tparam ibf_t Interleaved Bloom Filter type.
 * @param arguments Command line arguments.
 * @param ibf Interleaved Bloom Filter of the reference database.
 * @param thresholder Threshold for number of shared k-mers.
 * @param queue Shopping cart queue for load balancing between Valik prefiltering and Stellar search.
 */
template <typename ibf_t, typename TSequence>
void iterate_short_queries(search_arguments const & arguments,
                           ibf_t const & ibf,
                           raptor::threshold::threshold const & thresholder,
                           cart_queue<shared_query_record<TSequence>> & queue)
{
    using TId = seqan2::CharString;
    std::vector<shared_query_record<TSequence>> query_records{};
    constexpr uint64_t chunk_size = (1ULL << 20) * 10;

    seqan2::SeqFileIn inSeqs;
    if (!open(inSeqs, arguments.query_file.c_str()))
    {
        throw std::runtime_error("Failed to open " + arguments.query_file.string() + " file.");
    }

    std::set<TId> uniqueIds; // set of short IDs (cut at first whitespace)
    bool idsUnique = true;

    size_t seqCount{0};
    for (; !atEnd(inSeqs); ++seqCount)
    {
        TSequence seq{};
        TId id{};
        readRecord(id, seq, inSeqs);
        idsUnique &= stellar::_checkUniqueId(uniqueIds, id);

        query_records.emplace_back(std::move(seq), seqan2::toCString(std::move(id)));

        if (query_records.size() > chunk_size)
        {
            prefilter_queries_parallel<shared_query_record<TSequence>>(ibf, arguments, query_records, thresholder, queue);
            query_records.clear();
        }
    }

    if (!idsUnique)
        std::cerr << "WARNING: Non-unique query ids. Output can be ambiguous.\n";

    prefilter_queries_parallel<shared_query_record<TSequence>>(ibf, arguments, query_records, thresholder, queue);
}

/**
 * @brief Function that creates split query records from fasta file input and sends them for prefiltering.
 *
 * @tparam ibf_t Interleaved Bloom Filter type.
 * @param arguments Command line arguments.
 * @param ibf Interleaved Bloom Filter of the reference database.
 * @param thresholder Threshold for number of shared k-mers.
 * @param queue Shopping cart queue for load balancing between Valik prefiltering and Stellar search.
 * @param meta Metadata table for split query segments.
 */
template <typename ibf_t, typename TSequence>
void iterate_split_queries(search_arguments const & arguments,
                           ibf_t const & ibf,
                           raptor::threshold::threshold const & thresholder,
                           cart_queue<shared_query_record<TSequence>> & queue,
                           metadata & meta)
{
    using TId = seqan2::CharString;
    std::vector<shared_query_record<TSequence>> query_records{};
    constexpr uint64_t chunk_size = (1ULL << 20) * 10;

    seqan2::SeqFileIn inSeqs;
    if (!open(inSeqs, arguments.query_file.c_str()))
    {
        throw std::runtime_error("Failed to open " + arguments.query_file.string() + " file.");
    }

    std::set<TId> uniqueIds; // set of short IDs (cut at first whitespace)
    bool idsUnique = true;

    size_t seqCount{0};
    for (; !atEnd(inSeqs); ++seqCount)
    {
        TSequence seq{};
        TId id{};
        readRecord(id, seq, inSeqs);
        idsUnique &= stellar::_checkUniqueId(uniqueIds, id);

        auto query_ptr = std::make_shared<TSequence>(std::move(seq));

        for (auto const & seg : meta.segments_from_ind(seqCount))
        {
            // each split query record contains a copy of the same shared pointer
            query_records.emplace_back(seqan2::toCString(std::move(id)), seg, query_ptr);

            if (query_records.size() > chunk_size)
            {
                prefilter_queries_parallel<shared_query_record<TSequence>>(ibf, arguments, query_records, thresholder, queue);
                query_records.clear();  // shared pointers are erased -> memory is deallocated
            }
        }
    }

    if (!idsUnique)
        std::cerr << "WARNING: Non-unique query ids. Output can be ambiguous.\n";

    prefilter_queries_parallel<shared_query_record<TSequence>>(ibf, arguments, query_records, thresholder, queue);
}

}   // namespace valik::app
