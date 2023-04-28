#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <iomanip>

// Read two gene sequences from a file
void read_gene_sequences(const std::string &file_name, std::string &seq1, std::string &seq2) {
    std::ifstream file(file_name);
    if (file.is_open()) {
        std::getline(file, seq1);
        std::getline(file, seq2);
        file.close();
    }
}

// Compute the longest common subsequence matrix
std::vector<std::vector<int>> compute_lcs_matrix(const std::string &seq1, const std::string &seq2) {
    int m = seq1.length();
    int n = seq2.length();
    std::vector<std::vector<int>> c(m + 1, std::vector<int>(n + 1));

    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            if (seq1[i - 1] == seq2[j - 1]) {
                c[i][j] = c[i - 1][j - 1] + 1;
            } else {
                c[i][j] = std::max(c[i - 1][j], c[i][j - 1]);
            }
        }
    }
    return c;
}

// Reconstruct the longest common subsequence from the matrix
std::string reconstruct_lcs(const std::string &seq1, const std::string &seq2, const std::vector<std::vector<int>> &c) {
    int i = seq1.length();
    int j = seq2.length();
    std::string lcs;

    while (i > 0 && j > 0) {
        if (seq1[i - 1] == seq2[j - 1]) {
            lcs.push_back(seq1[i - 1]);
            i--;
            j--;
        } else if (c[i - 1][j] > c[i][j - 1]) {
            i--;
        } else {
            j--;
        }
    }
    std::reverse(lcs.begin(), lcs.end());
    return lcs;
}

// Read multiple gene sequences from a file
void read_multi_gene_sequences(const std::string &file_name, std::vector<std::string> &sequences) {
    std::ifstream file(file_name);
    if (file.is_open()) {
        int num_sequences;
        file >> num_sequences;
        file.ignore();

        for (int i = 0; i < num_sequences; ++i) {
            std::string seq;
            std::getline(file, seq);
            sequences.push_back(seq);
        }
        file.close();
    }
}

// Compute the length of the longest common subsequence using dynamic programming
int compute_lcs_length(const std::string &seq1, const std::string &seq2) {
    std::vector<std::vector<int>> c = compute_lcs_matrix(seq1, seq2);
    return c[seq1.length()][seq2.length()];
}

// Determine the similarity between two sequences based on their LCS length and sequence lengths
char determine_similarity(const std::string &seq1, const std::string &seq2) {
    int m = seq1.length();
    int n = seq2.length();
    int lcs_length = compute_lcs_length(seq1, seq2);
    int shorter_length = std::min(m, n);

    double length_ratio = static_cast<double>(shorter_length) / std::max(m, n);
    double lcs_ratio = static_cast<double>(lcs_length) / shorter_length;
    if (length_ratio >= 0.9 && lcs_ratio >= 0.9) {
        return 'H';
    } else if (length_ratio >= 0.8 && lcs_ratio >= 0.8) {
        return 'M';
    } else if (length_ratio >= 0.6 && lcs_ratio >= 0.5) {
        return 'L';
    } else {
        return 'D';
    }
}
// Output a table showing the similarity between each pair of sequences in the input vector
void output_similarity_table(const std::vector<std::string> &sequences) {
    int num_sequences = sequences.size();
    std::cout << std::setw(4);
    for (int i = 1; i <= num_sequences; i++) {
        std::cout << std::setw(3) << i;
    }
    std::cout << std::endl;

    for (int i = 0; i < num_sequences; i++) {
        std::cout << std::setw(2) << i + 1 << " ";
        for (int j = 0; j < i; j++) {
            std::cout << "   ";
        }
        for (int j = i; j < num_sequences; j++) {
            if (i == j) {
                std::cout << " - ";
            } else {
                char sim = determine_similarity(sequences[i], sequences[j]);
                std::cout << std::setw(3) << sim;
            }
        }
        std::cout << std::endl;
    }
}
int main() {
    std::string seq1, seq2;
    read_gene_sequences("twoStrings.txt", seq1, seq2);
    std::vector<std::vector<int>> c = compute_lcs_matrix(seq1, seq2);
    std::string lcs = reconstruct_lcs(seq1, seq2, c);
    std::cout << "Sequence 1: " << seq1 << std::endl;
    std::cout << "Sequence 2: " << seq2 << std::endl;
    std::cout << "Longest common subsequence: " << lcs << std::endl;
    std::vector<std::string> sequences;
    read_multi_gene_sequences("multiStrings.txt", sequences);

    output_similarity_table(sequences);

    return 0;
}