#include "./databases/FHE_disk_database.hpp"
#include "globals.hpp"
#include "./databases/tools.hpp"

#include <iostream>
#include <string>
#include <vector>

struct Config {
    int row = -1;
    int columns = -1;
    int threads = 1;
    std::string vcf_file;
    std::string csv_file;
};

int SEED = 27;

void print_help() {
    std::cout << "Usage:\n";
    std::cout << "  --rows <int>         Specify rows count\n";
    std::cout << "  --columns <int>     Specify columns count\n";
    std::cout << "  --threads <int>     Specify number of threads (default: 1)\n";
    std::cout << "  --vcf <filepath>    Specify path to a VCF file\n";
    std::cout << "  --csv <filepath>    Specify path to a CSV file\n";
    std::cout << "  --help              Display this help message\n";
    std::cout << "\nRequirements:\n";
    std::cout << "  - Either both --vcf and --csv must be provided, or both --rows and --columns.\n";
    std::cout << "  - You cannot mix row/columns with VCF/CSV options.\n";
}

bool parse_arguments(int argc, char* argv[], Config& config) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "--help") {
            print_help();
            return false;
        } else if (arg == "--rows" && i + 1 < argc) {
            config.row = std::stoi(argv[++i]);
        } else if (arg == "--columns" && i + 1 < argc) {
            config.columns = std::stoi(argv[++i]);
        } else if (arg == "--threads" && i + 1 < argc) {
            config.threads = std::stoi(argv[++i]);
        } else if (arg == "--vcf" && i + 1 < argc) {
            config.vcf_file = argv[++i];
        } else if (arg == "--csv" && i + 1 < argc) {
            config.csv_file = argv[++i];
        } else {
            std::cerr << "Unknown argument: " << arg << std::endl;
            print_help();
            return false;
        }
    }

    // Validation of requirements
    bool hasRowCol = (config.row != -1 && config.columns != -1);
    bool hasFiles = (!config.vcf_file.empty() && !config.csv_file.empty());

    if (!(hasRowCol || hasFiles) || (hasRowCol && hasFiles)) {
        std::cerr << "Error: Either both --vcf and --csv must be provided to run on real data, or both --row and --columns to generate synthetic data, but not both." << std::endl;
        print_help();
        return false;
    }

    return true;
}

int main(int argc, char* argv[]) {
    Config config;

    if (!parse_arguments(argc, argv, config)) {
        return 1;
    }

    // Output parsed values for debugging purposes
    std::cout << "Parsed arguments:\n";
    if (config.row != -1) std::cout << "  Row: " << config.row << "\n";
    if (config.columns != -1) std::cout << "  Columns: " << config.columns << "\n";
    std::cout << "  Threads: " << config.threads << "\n";
    if (!config.vcf_file.empty()) std::cout << "  VCF File: " << config.vcf_file << "\n";
    if (!config.csv_file.empty()) std::cout << "  CSV File: " << config.csv_file << "\n";

    // Initialize the database
    FHEDiskDatabase *dbFHEInstance = new FHEDiskDatabase(constants::Large, false);
    dbFHEInstance->setNumThreadsEncryption(config.threads);
    if (config.row != -1) {
        // Generate synthetic data
        vector<vector<int32_t>> db;

        for (int c = 0; c < config.columns; ++c) {
            vector<int32_t> column;
            for (int r = 0; r < config.row; ++r) {
                column.push_back(rand() % 3);
            }
            db.push_back(column);
        }
        dbFHEInstance->genData(config.row, config.columns, SEED);
    } else {
        // Load data from files
        dbFHEInstance->setVCFData(config.vcf_file);
        dbFHEInstance->setPhenoData(config.csv_file);     
    }

    std::cout << "Data loaded / generated successfully." << std::endl;

    return 0;
}