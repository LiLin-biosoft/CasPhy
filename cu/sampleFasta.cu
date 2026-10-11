#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <cstdlib>

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <input.fasta> <indices.txt> <output.fasta>\n";
        return 1;
    }

    const char* fasta_path  = argv[1];
    const char* index_path  = argv[2];
    const char* output_path = argv[3];

    // Read list of target indices from the index file
    std::ifstream ifs_idx(index_path);
    if (!ifs_idx) {
        std::cerr << "Error: cannot open index file '" << index_path << "'\n";
        return 1;
    }
    std::vector<int> targets;
    while (true) {
        std::string token;
        if (!(ifs_idx >> token)) break;
        try {
            int x = std::stoi(token);
            if (x > 0) {
                targets.push_back(x);
            } else {
                std::cerr << "Warning: skipping non-positive index " << x << "\n";
            }
        } catch (...) {
            std::cerr << "Warning: invalid index token '" << token << "' skipped\n";
        }
    }
    ifs_idx.close();

    if (targets.empty()) {
        std::cerr << "Error: no valid indices found in '" << index_path << "'\n";
        return 1;
    }

    // Sort, deduplicate, and build lookup set
    std::sort(targets.begin(), targets.end());
    targets.erase(std::unique(targets.begin(), targets.end()), targets.end());
    std::set<int> target_set(targets.begin(), targets.end());

    // Map to record which indices have been found
    std::map<int,bool> found;
    for (int idx : targets) {
        found[idx] = false;
    }

    // Open input FASTA and output file
    std::ifstream fin(fasta_path);
    if (!fin) {
        std::cerr << "Error: cannot open input FASTA '" << fasta_path << "'\n";
        return 1;
    }
    std::ofstream fout(output_path);
    if (!fout) {
        std::cerr << "Error: cannot open output file '" << output_path << "'\n";
        return 1;
    }

    std::string line;
    int current_index = 0;
    bool write_mode = false;
    int extracted = 0;

    // Read FASTA line by line
    while (std::getline(fin, line)) {
        if (!line.empty() && line[0] == '>') {
            // New sequence header encountered
            ++current_index;
            if (target_set.count(current_index)) {
                write_mode = true;
                found[current_index] = true;
                ++extracted;
                fout << line << "\n";
            } else {
                write_mode = false;
            }
        } else {
            // Sequence line
            if (write_mode) {
                fout << line << "\n";
            }
        }
    }

    fin.close();
    fout.close();

    // Report any indices that were not found
    for (int idx : targets) {
        if (!found[idx]) {
            std::cerr << "Warning: sequence #" << idx
                      << " not found in '" << fasta_path << "'\n";
        }
    }

    if (extracted > 0) {
        std::cout << "Successfully extracted " << extracted
                  << " sequences to '" << output_path << "'\n";
        return 0;
    } else {
        std::cerr << "Error: no sequences were extracted.\n";
        return 1;
    }
}
