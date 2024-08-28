#include "hf.h"

int main(int argc, char** argv) {
    auto[opt, timers] = Acorn::HF::initialize(argc, argv); Acorn::HF::run(opt, timers);
}
