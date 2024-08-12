#pragma once

#include "tensor.h"

class System {
public:
    System(const std::string& path); double nuclearRepulsion() const; int nocc() const;

private:
    std::vector<int> AN; torch::Tensor R;
};
