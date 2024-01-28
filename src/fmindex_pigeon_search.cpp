#include <sstream>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>


int main(int argc, char const* const* argv) {
    using namespace std;
    seqan3::argument_parser parser{"fmindex_pigeon_search", argc, argv, seqan3::update_notifications::off};

    parser.info.author = "SeqAn-Team";
    parser.info.version = "1.0.0";

    auto index_path = std::filesystem::path{};
    parser.add_option(index_path, '\0', "index", "path to the query file");

    auto reference_file = std::filesystem::path{};
    parser.add_option(reference_file, '\0', "reference", "path to the reference file");

    auto query_file = std::filesystem::path{};
    parser.add_option(query_file, '\0', "query", "path to the query file");

    auto number_of_queries = size_t{100};
    parser.add_option(number_of_queries, '\0', "query_ct", "number of query, if not enough queries, these will be duplicated");

    auto number_of_errors = uint8_t{0};
    parser.add_option(number_of_errors, '\0', "errors", "number of allowed hamming distance errors");

    try {
         parser.parse();
    } catch (seqan3::argument_parser_error const& ext) {
        seqan3::debug_stream << "Parsing error. " << ext.what() << "\n";
        return EXIT_FAILURE;
    }

    // loading our files
    auto reference_stream = seqan3::sequence_file_input{reference_file};
    auto query_stream     = seqan3::sequence_file_input{query_file};

    // read reference into memory
    std::vector<seqan3::dna5> reference;
    for (auto& record : reference_stream) {
        auto r = record.sequence();
        reference.insert(reference.end(), r.begin(), r.end());
    }

    // read query into memory
    std::vector<std::vector<seqan3::dna5>> queries;
    for (auto& record : query_stream) {
        queries.push_back(record.sequence());
    }

    // loading fm-index into memory
    using Index = decltype(seqan3::fm_index{std::vector<std::vector<seqan3::dna5>>{}}); // Some hack
    Index index; // construct fm-index
    {
        seqan3::debug_stream << "Loading 2FM-Index ... " << std::flush;
        std::ifstream is{index_path, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(index);
        seqan3::debug_stream << "done\n";
    }

    // duplicate input until its large enough
    while (queries.size() < number_of_queries) {
        auto old_count = queries.size();
        queries.resize(2 * old_count);
        std::copy_n(queries.begin(), old_count, queries.begin() + old_count);
    }
    queries.resize(number_of_queries); // will reduce the amount of searches

    using std::chrono::high_resolution_clock;
    auto start = high_resolution_clock::now();
    seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{0}};
    //!TODO !ImplementMe use the seqan3::search to find a partial error free hit, verify the rest inside the text

    for (auto & query : queries)
    {
        std::vector<std::vector<seqan3::dna5>> parts(number_of_errors+1);
        int segment_size = (query.size() / (number_of_errors+1)) + 1;

        for (int i = 0; i < query.size()-1; i++)
        {
            parts[i / segment_size].push_back(query[i]);
        }
        auto results = seqan3::search(parts, index, cfg);
        
        std::unordered_set<int> occurenceIndices;
        for (auto & result : results)
        {  
            int start = result.reference_begin_position() - (result.query_id()*segment_size);
            int end = start + query.size() - 1;
            if (end < reference.size() && start >= 0) {
                if (occurenceIndices.count(start) == 0) {
                    std::vector<seqan3::dna5> reference_segment;
                    for (int i = start; i < end; i++) {
                        reference_segment.push_back(reference[i]);
                    }

                    int distance = 0;
                    for (int i = 1; i < reference_segment.size(); i++)
                    {
                        if (reference_segment[i] != query[i]) {
                            distance++;
                        }
                    }

                    if (distance < number_of_errors)
                    {
                        occurenceIndices.insert(start);
                    }
                }
            }
        }
    }
    auto end = high_resolution_clock::now();  
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    std::cout << "Time took: " << duration.count() << " ns\n";
    return 0;
}
