.RECIPEPREFIX = >

ARCH := x86_64
OS   :=  linux
ABI  :=   musl

ZIG_VERSION      := 0.14.1
ZLS_VERSION      := 0.14.0
EIGEN_VERSION    :=  3.4.0
LIBINT_VERSION   := 2.11.1
LLVM_VERSION     := 20.1.8
OPENBLAS_VERSION := 0.3.30

EXPRTK_COMMIT := cc1b800c2bd1ac3ac260478c915d2aec6f4eb41c

URL_EIGEN    :=                             https://gitlab.com/libeigen/eigen/-/archive/$(EIGEN_VERSION)/eigen-$(EIGEN_VERSION).tar.gz
URL_LIBINT   :=                    https://github.com/evaleev/libint/releases/download/v$(LIBINT_VERSION)/libint-$(LIBINT_VERSION).tgz
URL_LLVM     := https://github.com/llvm/llvm-project/releases/download/llvmorg-$(LLVM_VERSION)/llvm-project-$(LLVM_VERSION).src.tar.xz
URL_OPENBLAS :=     https://github.com/OpenMathLib/OpenBLAS/releases/download/v$(OPENBLAS_VERSION)/OpenBLAS-$(OPENBLAS_VERSION).tar.gz
URL_EXPRTK   :=                                       https://raw.githubusercontent.com/ArashPartow/exprtk/$(EXPRTK_COMMIT)/exprtk.hpp
URL_ZIG      :=                                    https://ziglang.org/download/$(ZIG_VERSION)/zig-$(ARCH)-$(OS)-$(ZIG_VERSION).tar.xz
URL_ZLS      :=                              https://github.com/zigtools/zls/releases/download/$(ZLS_VERSION)/zls-$(ARCH)-$(OS).tar.xz

export OMP_NUM_THREADS := 1

all: acorn

# ACORN BUILDING TARGETS ===============================================================================================================================================================================

.PHONY: acorn docs test

acorn: zig-bin/.done external-$(ARCH)-$(OS)/.done
> ./zig-bin/zig build --summary all

docs: zig-bin/.done external-$(ARCH)-$(OS)/.done
> ./zig-bin/zig build --summary all docs

test: zig-bin/.done external-$(ARCH)-$(OS)/.done
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
> @mkdir -p external-$(ARCH)-$(OS)/include external-$(ARCH)-$(OS)/lib zig-bin

# LIBRARY TARGETS ======================================================================================================================================================================================

LIBRARIES := eigen.tar.gz libint.tar.gz llvm.tar.xz openblas.tar.gz include/exprtk.hpp

external-$(ARCH)-$(OS)/.done: $(addprefix external-$(ARCH)-$(OS)/,$(LIBRARIES))
> ./script/library.sh && touch $@

external-$(ARCH)-$(OS)/eigen.tar.gz:       URL =    $(URL_EIGEN)
external-$(ARCH)-$(OS)/libint.tar.gz:      URL =   $(URL_LIBINT)
external-$(ARCH)-$(OS)/llvm.tar.xz:        URL =     $(URL_LLVM)
external-$(ARCH)-$(OS)/openblas.tar.gz:    URL = $(URL_OPENBLAS)
external-$(ARCH)-$(OS)/include/exprtk.hpp: URL =   $(URL_EXPRTK)

external-$(ARCH)-$(OS)/%: | tree
> wget -q -O $@ $(URL)

# COMPILER TARGETS =====================================================================================================================================================================================

zig-bin/.done: zig-bin/zig zig-bin/zls
> touch $@

zig-bin/zig: | tree
> wget -q -O zig.tar.xz $(URL_ZIG) && tar -xf zig.tar.xz -C zig-bin --strip-components=1 && rm zig.tar.xz

zig-bin/zls: | tree
> wget -q -O zls.tar.xz $(URL_ZLS) && tar -xf zls.tar.xz -C zig-bin && rm zls.tar.xz

# CLEAN TARGETS ========================================================================================================================================================================================

clean:
> git clean -dffx

clean-cache:
> rm -rf .zig-cache

clean-output:
> rm -rf zig-out
