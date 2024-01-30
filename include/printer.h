#pragma once

#include "eigen.h"

struct Printer {
    static void Title(const std::string& title, bool newline = true);
    static void Print(const Tensor<>& A, const std::string& title);
    static void Print(const Matrix<>& A, const std::string& title);
    static void Print(double a, const std::string& title);
};

#include <iostream>
