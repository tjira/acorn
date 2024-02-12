#include "printer.h"

void Printer::Print(double a, const std::string& title) {
    std::printf("%s: %.14f\n", title.c_str(), a);
}

void Printer::Print(const Matrix<>& A, const std::string& title) {
    std::cout << title << ":\n" << A; std::printf("\n%s NORM AND SUM: %.2e %.2e\n", title.c_str(), A.norm(), A.sum());
}

void Printer::Print(const Tensor<3>& A, const std::string& title) {
    Print(Eigen::Map<const Matrix<>>(A.data(), A.dimension(0) * A.dimension(2), A.dimension(1)), title);
}

void Printer::Print(const Tensor<4>& A, const std::string& title) {
    Print(Eigen::Map<const Matrix<>>(A.data(), A.dimension(0) * A.dimension(2), A.dimension(1) * A.dimension(3)), title);
}

void Printer::Print(const Tensor<5>& A, const std::string& title) {
    Print(Eigen::Map<const Matrix<>>(A.data(), A.dimension(0) * A.dimension(2) * A.dimension(4), A.dimension(1) * A.dimension(3)), title);
}

void Printer::Title(const std::string& title, bool newline) {
    std::cout << std::string(146, '-') + "\n" << title << "\n" << std::string(146, '-') << std::endl; if (newline) std::cout << std::endl;
}
