#include "expression.h"

int main(int argc, char** argv) {
    auto [opt, timers] = Acorn::Expression::initialize(argc, argv); Acorn::Expression::run(opt, timers);
}
