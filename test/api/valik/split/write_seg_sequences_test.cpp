#include <gtest/gtest.h>

#include <valik/split/write_seg_sequences.hpp>

#include <seqan3/alphabet/nucleotide/dna4.hpp>

// Generate the full path of a test input file that is provided in the data directory.
static std::filesystem::path data(std::string const & filename)
{
    return std::filesystem::path{std::string{DATADIR}}.concat(filename);
}

static inline std::string const string_from_file(std::filesystem::path const & path, std::ios_base::openmode const mode = std::ios_base::in)
{
    std::ifstream file_stream(path, mode);
    if (!file_stream.is_open())
        throw std::logic_error{"Cannot open " + path.string()};
    std::stringstream file_buffer;
    file_buffer << file_stream.rdbuf();
    return {file_buffer.str()};
}

static void trim_fasta_id(std::string & seg_seq)
{
    size_t pos = seg_seq.find('\n');
    seg_seq = seg_seq.substr(pos);
    seg_seq.erase(std::remove_if(seg_seq.begin(),
                seg_seq.end(),
                [](unsigned char x){ return std::isspace(x); }),
                seg_seq.end());
}

static void const test_write_out(size_t overlap, size_t bins)
{
    std::string path_prefix = "write_out_" + std::to_string(overlap) + "_" + std::to_string(bins);

    valik::reference_metadata reference(data(path_prefix + "_reference_metadata.txt"), false);
    valik::reference_segments segments(data(path_prefix + "_reference_segments.txt"));
    valik::write_seg_sequences<seqan3::dna4>(reference, segments, data(path_prefix + "_ref.fasta"));

    for (size_t i = 0; i < bins - 1; i++)
    {
        valik::reference_segments::segment current_seg = segments.members[i];
        valik::reference_segments::segment next_seg = segments.members[i + 1];

        std::string current_seg_seq = string_from_file(data(path_prefix + "_ref_" + std::to_string(i) + ".fasta"), std::ios::binary);
        std::string next_seg_seq = string_from_file(data(path_prefix + "_ref_" + std::to_string(i + 1) + ".fasta"), std::ios::binary);
        trim_fasta_id(current_seg_seq);
        trim_fasta_id(next_seg_seq);

        EXPECT_EQ(current_seg_seq.size(), current_seg.len);
        EXPECT_EQ(next_seg_seq.size(), next_seg.len);

        if (current_seg.ref_ind == next_seg.ref_ind)
        {
            EXPECT_EQ(current_seg_seq.substr(current_seg_seq.size() - overlap), next_seg_seq.substr(0, overlap));
        }
    }
}

TEST(write_seg_sequences, o0_b4)
{
    size_t overlap = 0;
    size_t bins = 4;

    test_write_out(overlap, bins);
}

TEST(write_seg_sequences, o20_b4)
{
    size_t overlap = 20;
    size_t bins = 4;

    test_write_out(overlap, bins);
}

TEST(write_seg_sequences, o0_b16)
{
    size_t overlap = 0;
    size_t bins = 16;

    test_write_out(overlap, bins);
}

TEST(write_seg_sequences, o20_b16)
{
    size_t overlap = 20;
    size_t bins = 16;

    test_write_out(overlap, bins);
}
