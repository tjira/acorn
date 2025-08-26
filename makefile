.RECIPEPREFIX = >

DEBUG ?= 0
NPROC ?= 1

ARCH := $(shell uname -m | tr '[:upper:]' '[:lower:]')
OS   := $(shell uname -s | tr '[:upper:]' '[:lower:]')

TARGETS := aarch64-linux x86_64-linux

ZIG_VERSION      := 0.15.1
ZLS_VERSION      := 0.15.0
EIGEN_VERSION    :=  3.4.0
LIBINT_VERSION   := 2.11.1
LLVM_VERSION     := 20.1.8
OPENBLAS_VERSION := 0.3.30

EXPRTK_COMMIT := cc1b800c2bd1ac3ac260478c915d2aec6f4eb41c

export OMP_NUM_THREADS := 1

all: acorn

# ACORN BUILDING TARGETS ===============================================================================================================================================================================

.PHONY: acorn docs test

acorn: $(addprefix external-,$(addsuffix /.done,$(TARGETS)))
> ./zig-bin/zig build $(if $(filter 1,$(DEBUG)),-DDEBUG) -j$(NPROC)

docs: $(addprefix external-,$(addsuffix /.done,$(TARGETS)))
> ./zig-bin/zig build $(if $(filter 1,$(DEBUG)),-DDEBUG) -j$(NPROC) docs

run: $(addprefix external-,$(addsuffix /.done,$(TARGETS)))
> ./zig-bin/zig build $(if $(filter 1,$(DEBUG)),-DDEBUG) -j$(NPROC) run

test: $(addprefix external-,$(addsuffix /.done,$(TARGETS)))
> ./zig-bin/zig build $(if $(filter 1,$(DEBUG)),-DDEBUG) -j$(NPROC) test

# BENCHMARKING TARGETS =================================================================================================================================================================================

BENCHMARK_TARGETS := contract dgees dgesdd dgemm dsyevd

BENCHMARK_ARGS.contract := contract  100
BENCHMARK_ARGS.dgees    := dgees    1000
BENCHMARK_ARGS.dgesdd   := dgesdd   1000
BENCHMARK_ARGS.dgemm    := dgemm    3000
BENCHMARK_ARGS.dsyevd   := dsyevd   2000

BENCHMARK_SAMPLES := 1

benchmark: $(addprefix benchmark-,$(BENCHMARK_TARGETS))

benchmark-%: acorn
> ./zig-out/$(ARCH)-$(OS)-musl/benchmark $(BENCHMARK_ARGS.$*) $(BENCHMARK_SAMPLES)

# CUSTOM TARGETS --=====================================================================================================================================================================================

lines:
> @git ls-files | grep -E '\.cpp$$|\.h$$|\.zig$$' | xargs cat | awk 'NF' | wc -l

linguist:
> @github-linguist

serve:
> @cd docs && bundle exec jekyll serve

# EXTERNAL TARGETS =====================================================================================================================================================================================

EXTERNAL := eigen.tar.gz libint.tar.gz llvm.tar.xz openblas.tar.gz include/exprtk.hpp

$(addprefix external-,$(addsuffix /.done,$(TARGETS))): external-%/.done: $(addprefix external-$(ARCH)-$(OS)/,$(EXTERNAL)) zig-bin/zig
> ./script/library.sh $* && touch $@

external-$(ARCH)-$(OS)/eigen.tar.gz:       URL =                             https://gitlab.com/libeigen/eigen/-/archive/$(EIGEN_VERSION)/eigen-$(EIGEN_VERSION).tar.gz
external-$(ARCH)-$(OS)/libint.tar.gz:      URL =                    https://github.com/evaleev/libint/releases/download/v$(LIBINT_VERSION)/libint-$(LIBINT_VERSION).tgz
external-$(ARCH)-$(OS)/llvm.tar.xz:        URL = https://github.com/llvm/llvm-project/releases/download/llvmorg-$(LLVM_VERSION)/llvm-project-$(LLVM_VERSION).src.tar.xz
external-$(ARCH)-$(OS)/openblas.tar.gz:    URL =     https://github.com/OpenMathLib/OpenBLAS/releases/download/v$(OPENBLAS_VERSION)/OpenBLAS-$(OPENBLAS_VERSION).tar.gz
external-$(ARCH)-$(OS)/include/exprtk.hpp: URL =                                       https://raw.githubusercontent.com/ArashPartow/exprtk/$(EXPRTK_COMMIT)/exprtk.hpp

external-$(ARCH)-$(OS)/%:
> mkdir -p external-$(ARCH)-$(OS)/include && wget -q -O $@ $(URL)

zig-bin/zig:
> mkdir -p zig-bin && wget -q -O - https://ziglang.org/download/$(ZIG_VERSION)/zig-$(ARCH)-$(OS)-$(ZIG_VERSION).tar.xz       | tar -Jx -C zig-bin --strip-components=1 && touch zig-bin/zig
> mkdir -p zig-bin && wget -q -O - https://github.com/zigtools/zls/releases/download/$(ZLS_VERSION)/zls-$(ARCH)-$(OS).tar.xz | tar -Jx -C zig-bin --strip-components=0 && touch zig-bin/zls

# CLEAN TARGETS ========================================================================================================================================================================================

clean:
> git clean -dffx

clean-binary:
> rm -rf zig-out

clean-cache:
> rm -rf ${HOME}/.cache/zig .zig-cache

clean-root:
> rm -rf *.all *.json *.mat *.xyz .bagel .orca
