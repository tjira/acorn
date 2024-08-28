#include "cc.h"

int main(int argc, char** argv) {
    auto[opt, timers] = Acorn::CC::initialize(argc, argv); Acorn::CC::run(opt, timers);
}
