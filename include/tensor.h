#pragma once

#include <torch/torch.h>
#include <fstream>

namespace torch {
    template <int D> torch::Tensor ReadTensor(const std::string& path); template <int D> void WriteTensor(const std::string& path, const torch::Tensor& ten);
}

template <int D>
torch::Tensor torch::ReadTensor(const std::string& path) {
    // open the input file and check for errors
    std::ifstream file(path); if (!file.good()) throw std::runtime_error("UNABLE TO OPEN FILE `" + path + "` FOR READING");

    // read the dimensions and assign them to an array
    std::vector<int64_t> dims(D); size_t prod = 1; for (int i = 0; i < D; i++) file >> dims.at(i), prod *= dims.at(i);

    // define the data container and read the file
    std::vector<double> data(prod); for (size_t i = 0; i < prod; i++) file >> data.at(i);

    // return the tensor
    return torch::from_blob(data.data(), dims, torch::TensorOptions().dtype(torch::kDouble)).clone();
}

template <int D>
void torch::WriteTensor(const std::string& path, const torch::Tensor& ten) {
    // open the input file and check for errors
    std::ofstream file(path); if (!file.good()) throw std::runtime_error("UNABLE TO OPEN FILE `" + path + "` FOR WRITING");

    // write the dimensions to the header and set the formatting
    for (size_t i = 0; i < ten.sizes().size(); i++) {file << ten.sizes().at(i) << (i < ten.sizes().size() - 1 ? " " : "");} file << "\n" << std::fixed << std::setprecision(14);

    // write 2nd order tensor
    if constexpr (D == 2) {
        for (int i = 0; i < ten.sizes().at(0); i++) {
            for (int j = 0; j < ten.sizes().at(1); j++) {
                file << std::setw(20) << ten.index({i, j}).item().toDouble() << (j < ten.sizes().at(1) - 1 ? " " : "");
            } file << "\n";
        }
    }

    // write 4th order tensor
    if constexpr (D == 4) {
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
