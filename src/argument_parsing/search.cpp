#include <seqan3/io/views/async_input_buffer.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <sliding_window/argument_parsing/search.hpp>
#include <sliding_window/search/search.hpp>

namespace sliding_window::app
{

void init_search_parser(seqan3::argument_parser & parser, search_arguments & arguments)
{
    init_shared_meta(parser);
    init_shared_options(parser, arguments);
    parser.add_option(arguments.ibf_file,
                      '\0',
                      "index",
                      "Provide a valid path to an IBF.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{});
    parser.add_option(arguments.query_file,
                      '\0',
                      "query",
                      "Provide a path to the query file.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{});
    parser.add_option(arguments.out_file,
                      '\0',
                      "output",
                      "Please provide a valid path to the output.",
                      seqan3::option_spec::required,
                      seqan3::output_file_validator{});
    parser.add_option(arguments.errors,
                      '\0',
                      "error",
                      "Choose the number of errors.",
                      seqan3::option_spec::standard,
                      positive_integer_validator{true});
    parser.add_option(arguments.tau,
                      '\0',
                      "tau",
                      "Threshold for probabilistic models.",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{0, 1});
    parser.add_option(arguments.threshold,
                      '\0',
                      "threshold",
                      "If set, this threshold is used instead of the probabilistic models.",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{0, 1});
    parser.add_option(arguments.pattern_size,
                      '\0',
                      "pattern",
                      "Choose the pattern size. Default: Use median of sequence lengths in query file.",
                      seqan3::option_spec::standard);
    parser.add_option(arguments.overlap,
                      '\0',
                      "overlap",
                      "Choose how much sequential sliding windows overlap. Default: pattern size - 1.",
                      seqan3::option_spec::standard);
    parser.add_flag(arguments.compressed,
                    '\0',
                    "compressed",
                    "Build a compressed IBF.");
    parser.add_flag(arguments.write_time,
                    '\0',
                    "time",
                    "Write timing file.",
                    seqan3::option_spec::advanced);
}

// TODO: split creating seqan3::argument_parser and checks/validation of the arguments
void run_search(seqan3::argument_parser & parser)
{
    search_arguments arguments{};

    init_search_parser(parser, arguments);

    try_parsing(parser);

    // ==========================================
    // Various checks.
    // ==========================================

    seqan3::input_file_validator<seqan3::sequence_file_input<>>{}(arguments.query_file);
    arguments.treshold_was_set = parser.is_option_set("threshold");


    // ==========================================
    // Read window and kmer size.
    // ==========================================
    {
        std::ifstream is{arguments.ibf_file, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(arguments.kmer_size);
        iarchive(arguments.window_size);
    }

    // ==========================================
    // Process --pattern.
    // ==========================================
    
    if (parser.is_option_set("pattern"))
    {
        if (arguments.pattern_size < arguments.window_size)
            throw seqan3::argument_parser_error{"The minimiser window cannot be bigger than the sliding window."};
    }
    else
        
	// Default behaviour: local match
	//arguments.pattern_size = arguments.window_size * 2;
    
	// Default behaviour: semiglobal match 
	if (!arguments.pattern_size)
	{
            std::vector<uint64_t> sequence_lengths{};
            seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>> query_in{arguments.query_file};
            for (auto & [seq] : query_in | seqan3::views::async_input_buffer(16))
            {
                sequence_lengths.push_back(std::ranges::size(seq));
            }
            std::sort(sequence_lengths.begin(), sequence_lengths.end());
            arguments.pattern_size = sequence_lengths[sequence_lengths.size()/2];
        }

    // ==========================================
    // More checks.
    // ==========================================

    if (parser.is_option_set("overlap"))
    {
	if (arguments.overlap >= arguments.pattern_size)
            throw seqan3::argument_parser_error{"The overlap size has to be smaller"
		    "than the sliding window (i.e pattern) size."};
    }
    else
        arguments.overlap = arguments.pattern_size - 1;

    // ==========================================
    // Dispatch
    // ==========================================
    sliding_window_search(arguments);
};

} // namespace sliding_window::app
