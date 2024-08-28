#include "cdyn.h"

int main(int argc, char** argv) {
    auto [opt, timers] = Acorn::CDYN::initialize(argc, argv); Acorn::CDYN::run(opt, timers);
}
