#pragma once

#include <torch/torch.h>

#include <fstream>
#include <iomanip>

#define OCC "abcdefghijklmnopqrstuvwxyz"
#define VRT "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

#define FORMAT(T) [](long ms) {char s[99]; std::sprintf(s, "%02ld:%02ld:%02ld.%03ld", ms / 3600000, ms % 3600000 / 60000, ms % 60000 / 1000, ms % 1000); return std::string(s);}(T)

typedef std::chrono::time_point<std::chrono::high_resolution_clock> timepoint; using namespace std::chrono; using namespace torch::indexing;
inline std::string eltime(timepoint t) {return FORMAT(duration_cast<milliseconds>(std::chrono::high_resolution_clock().now() - t).count());}

namespace torch {
    torch::Tensor ReadTensor(const std::string& path); void WriteTensor(const std::string& path, const torch::Tensor& ten);
}

inline torch::Tensor torch::ReadTensor(const std::string& path) {
    // open the input file and check for errors
    std::ifstream file(path); if (!file.good()) throw std::runtime_error("UNABLE TO OPEN FILE `" + path + "` FOR READING");

    // define the dimensions, its product and read the first line to a container
    std::vector<int64_t> dims; size_t prod = 1; std::string c; std::getline(file, c);

    // read the dimensions and calculate its product
    std::istringstream iss(c); while (iss >> c) dims.push_back(std::stoi(c)), prod *= dims.back();

    // define the data container and read the file
    std::vector<double> data(prod); for (size_t i = 0; i < prod; i++) file >> data.at(i);

    // return the tensor
    return torch::from_blob(data.data(), dims, torch::kDouble).clone();
}

inline void torch::WriteTensor(const std::string& path, const torch::Tensor& ten) {
    // open the input file and check for errors
    std::ofstream file(path); if (!file.good()) throw std::runtime_error("UNABLE TO OPEN FILE `" + path + "` FOR WRITING");

    // write the dimensions to the header and set the formatting
    for (size_t i = 0; i < ten.sizes().size(); i++) {file << ten.sizes().at(i) << (i < ten.sizes().size() - 1 ? " " : "");} file << "\n" << std::fixed << std::setprecision(14);

    // write 1st order tensor
    if (ten.sizes().size() == 1) {
        for (int i = 0; i < ten.sizes().at(0); i++) {
            file << std::setw(20) << ten.index({i}).item().toDouble() << "\n";
        }
    }

    // write 2nd order tensor
    if (ten.sizes().size() == 2) {
        for (int i = 0; i < ten.sizes().at(0); i++) {
            for (int j = 0; j < ten.sizes().at(1); j++) {
                file << std::setw(20) << ten.index({i, j}).item().toDouble() << (j < ten.sizes().at(1) - 1 ? " " : "");
            } file << "\n";
        }
    }

    // write 4th order tensor
    if (ten.sizes().size() == 4) {
        for (int i = 0; i < ten.sizes().at(0); i++) {
            for (int j = 0; j < ten.sizes().at(1); j++) {
                for (int k = 0; k < ten.sizes().at(2); k++) {
                    for (int l = 0; l < ten.sizes().at(3); l++) {
                        file << std::setw(20) << ten.index({i, j, k, l}).item().toDouble() << (k < ten.sizes().at(2) - 1 ? " " : "");
                    }
                } file << "\n";
            }
        }
    }
}
