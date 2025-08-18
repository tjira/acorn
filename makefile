.RECIPEPREFIX = >

ARCH := x86_64
OS   :=  linux
ABI  :=   musl

EXTERNAL_DIR := external-$(ARCH)-$(OS)

ZIG_VERSION      := 0.14.1
ZLS_VERSION      := 0.14.0
EIGEN_VERSION    :=  3.4.0
LIBINT_VERSION   := 2.11.1
LLVM_VERSION     := 20.1.8
OPENBLAS_VERSION := 0.3.30

EXPRTK_COMMIT := cc1b800c2bd1ac3ac260478c915d2aec6f4eb41c

export OMP_NUM_THREADS := 1

all: acorn

# ACORN BUILDING TARGETS ===============================================================================================================================================================================

.PHONY: acorn docs test

acorn: $(EXTERNAL_DIR)/.done
> ./zig-bin/zig build --summary all

docs: $(EXTERNAL_DIR)/.done
> ./zig-bin/zig build --summary all docs

test: $(EXTERNAL_DIR)/.done
> ./zig-bin/zig build --summary all test

# BENCHMARKING TARGETS =================================================================================================================================================================================

BENCHMARK_TARGETS := contract dgees dgemm dsyevd

BENCHMARK_ARGS.contract := contract  100
BENCHMARK_ARGS.dgees    := dgees    1000
BENCHMARK_ARGS.dgemm    := dgemm    3000
BENCHMARK_ARGS.dsyevd   := dsyevd   2000

BENCHMARK_SAMPLES := 1

benchmark: $(addprefix benchmark-,$(BENCHMARK_TARGETS))

benchmark-%: acorn
> ./zig-out/$(ARCH)-$(OS)-$(ABI)/benchmark $(BENCHMARK_ARGS.$*) $(BENCHMARK_SAMPLES)

# CUSTOM TARGETS --=====================================================================================================================================================================================

lines:
> @git ls-files | grep -E '\.cpp$$|\.h$$|\.zig$$' | xargs cat | awk 'NF' | wc -l

linguist:
> @github-linguist

serve:
> @cd docs && bundle exec jekyll serve

tree:
> @mkdir -p $(EXTERNAL_DIR)/include $(EXTERNAL_DIR)/lib zig-bin

# EXTERNAL TARGETS =====================================================================================================================================================================================

EXTERNAL := eigen.tar.gz libint.tar.gz llvm.tar.xz openblas.tar.gz zig.tar.xz zls.tar.xz include/exprtk.hpp

$(EXTERNAL_DIR)/.done: $(addprefix $(EXTERNAL_DIR)/,$(EXTERNAL))
> ./script/library.sh && touch $@

$(EXTERNAL_DIR)/eigen.tar.gz:       URL =                             https://gitlab.com/libeigen/eigen/-/archive/$(EIGEN_VERSION)/eigen-$(EIGEN_VERSION).tar.gz
$(EXTERNAL_DIR)/libint.tar.gz:      URL =                    https://github.com/evaleev/libint/releases/download/v$(LIBINT_VERSION)/libint-$(LIBINT_VERSION).tgz
$(EXTERNAL_DIR)/llvm.tar.xz:        URL = https://github.com/llvm/llvm-project/releases/download/llvmorg-$(LLVM_VERSION)/llvm-project-$(LLVM_VERSION).src.tar.xz
$(EXTERNAL_DIR)/openblas.tar.gz:    URL =     https://github.com/OpenMathLib/OpenBLAS/releases/download/v$(OPENBLAS_VERSION)/OpenBLAS-$(OPENBLAS_VERSION).tar.gz
$(EXTERNAL_DIR)/zig.tar.xz:         URL =                                    https://ziglang.org/download/$(ZIG_VERSION)/zig-$(ARCH)-$(OS)-$(ZIG_VERSION).tar.xz
$(EXTERNAL_DIR)/zls.tar.xz:         URL =                              https://github.com/zigtools/zls/releases/download/$(ZLS_VERSION)/zls-$(ARCH)-$(OS).tar.xz
$(EXTERNAL_DIR)/include/exprtk.hpp: URL =                                       https://raw.githubusercontent.com/ArashPartow/exprtk/$(EXPRTK_COMMIT)/exprtk.hpp

$(EXTERNAL_DIR)/%: | tree
> wget -q -O $@ $(URL)

# CLEAN TARGETS ========================================================================================================================================================================================

clean:
> git clean -dffx

clean-cache:
> rm -rf .zig-cache

clean-output:
> rm -rf zig-out
