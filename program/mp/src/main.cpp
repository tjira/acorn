#include "mp.h"

int main(int argc, char** argv) {
    auto[opt, timers] = Acorn::MBPT::initialize(argc, argv); Acorn::MBPT::run(opt, timers);
}
