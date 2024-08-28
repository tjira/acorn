#include "integral.h"

int main(int argc, char** argv) {
    auto [opt, timers] = Acorn::Integral::initialize(argc, argv); Acorn::Integral::run(opt, timers);
}
