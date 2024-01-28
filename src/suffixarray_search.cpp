#include <divsufsort.h>
#include <sstream>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

int main(int argc, char const* const* argv) {
    seqan3::argument_parser parser{"suffixarray_search", argc, argv, seqan3::update_notifications::off};

    parser.info.author = "SeqAn-Team";
    parser.info.version = "1.0.0";

    auto reference_file = std::filesystem::path{};
    parser.add_option(reference_file, '\0', "reference", "path to the reference file");

    auto query_file = std::filesystem::path{};
    parser.add_option(query_file, '\0', "query", "path to the query file");

    auto number_of_queries = size_t{100};
    parser.add_option(number_of_queries, '\0', "query_ct", "number of query, if not enough queries, these will be duplicated");

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
    // Attention: we are concatenating all sequences into one big combined sequence
    //            this is done to simplify the implementation of suffix_arrays
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

    // duplicate input until its large enough
    while (queries.size() < number_of_queries) {
        auto old_count = queries.size();
        queries.resize(2 * old_count);
        std::copy_n(queries.begin(), old_count, queries.begin() + old_count);
    }
    queries.resize(number_of_queries); // will reduce the amount of searches

    // Array that should hold the future suffix array
    std::vector<saidx_t> suffixarray;
    suffixarray.resize(reference.size()); // resizing the array, so it can hold the complete SA

    //!TODO !ImplementMe implement suffix array sort
    //Hint, if can use libdivsufsort (already integrated in this repo)
    //      https://github.com/y-256/libdivsufsort
    //      To make the `reference` compatible with libdivsufsort you can simply
    //      cast it by calling:
    //      `sauchar_t const* str = reinterpret_cast<sauchar_t const*>(reference.data());`
    sauchar_t const* str = reinterpret_cast<sauchar_t const*>(reference.data());
    divsufsort(str, suffixarray.data(), reference.size());

    using namespace std;
    using std::chrono::high_resolution_clock;
    auto start = high_resolution_clock::now();
    
    for (auto& q : queries) {
        //!TODO !ImplementMe apply binary search and find q  in reference using binary search on `suffixarray`
        // You can choose if you want to use binary search based on "naive approach", "mlr-trick", "lcp"
        vector<int> occurenceIndices;

        int left = 0;
        int right = reference.size() - 1;

        while (left <= right) {
            int middle = (left + right) / 2;
            bool flag = true;
            for (int i = 0; i < q.size()-1; i++) {
                if (q[i] != reference[suffixarray[middle] + i]) {
                    flag = false;
                    break;
                }
            }
            if (flag) {
                occurenceIndices.push_back(left);
                int middle_left = middle - 1;
                while (middle_left >= 0) {
                    bool flag = true;
                    for (int i = 0; i < q.size(); ++i) {
                        if (q[i] != reference[suffixarray[middle_left] + i]) {
                            flag = false;
                            break;
                        }
                    }
                    if (flag) {
                        occurenceIndices.push_back(middle_left);
                        middle_left--;
                    } else {
                        break;
                    }
                }
                int middle_right = middle + 1;
                while (middle_right < reference.size()) {
                    bool flag = true;
                    for (int i = 0; i < q.size()-1; i++) {
                        if (q[i] != reference[suffixarray[middle_right] + i]) {
                            flag = false;
                            break;
                        }
                    }
                    if (flag) {
                        occurenceIndices.push_back(middle_right);
                        middle_right++;
                    } else {
                        break;
                    }
                }
                break;
            } else {
                bool flag = true;
                for (int i = 0; i < q.size()-1; i++) {
                    if (q[i] != reference[suffixarray[middle] + i]) {
                        if (q[i] > reference[suffixarray[middle] + i]) {
                            flag = false;
                        }
                        break;
                    }
                }
                if (flag) {
                    right = middle - 1;
                } else {
                    left = middle + 1;
                }
            }
        }
    }
    auto end = high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    std::cout << "Time Took: " << duration.count() << " ns\n";
    return 0;
}
