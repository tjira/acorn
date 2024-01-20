#include "printer.h"

void Printer::Print(double a, const std::string& title) {
    std::printf("%s: %.14f\n", title.c_str(), a);
}

void Printer::Print(const Matrix<>& A, const std::string& title) {
    std::cout << title << ":\n" << A; std::printf("\n%s NORM: %.2e\n", title.c_str(), A.norm());
}

void Printer::Print(const Tensor<>& A, const std::string& title) {
    Print(Eigen::Map<const Matrix<>>(A.data(), A.dimension(0) * A.dimension(2), A.dimension(1) * A.dimension(3)), title);
}
