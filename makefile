.RECIPEPREFIX = >

ZIG_VERSION = 0.14.1
ZLS_VERSION = 0.14.0

ZIG_FLAGS = --cache-dir zig-cache --global-cache-dir zig-global

EIGEN_VERSION    =  3.4.0
LIBINT_VERSION   = 2.11.1
LLVM_VERSION     = 20.1.8
OPENBLAS_VERSION = 0.3.30

BENCHMARK_SAMPLES = 1

all: acorn

# ACORN BUILDING TARGETS ===============================================================================================================================================================================

acorn: library
> ./zig-bin/zig build $(ZIG_FLAGS) -DDEBUG

full: library
> ./zig-bin/zig build $(ZIG_FLAGS) -DFULL

test: library
> ./zig-bin/zig build $(ZIG_FLAGS) test

# BENCHMARKING TARGETS =================================================================================================================================================================================

.PHONY: benchmark

benchmark: benchmark-contract benchmark-dgees benchmark-dgemm benchmark-dsyevd

benchmark-contract: full
> ./zig-out/x86_64-linux-musl/benchmark/contract 100 $(BENCHMARK_SAMPLES) | grep AVERAGE | tee -a benchmark.out

benchmark-dgees: full
> ./zig-out/x86_64-linux-musl/benchmark/dgees 1000 $(BENCHMARK_SAMPLES) | grep AVERAGE | tee -a benchmark.out

benchmark-dgemm: full
> ./zig-out/x86_64-linux-musl/benchmark/dgemm 3000 $(BENCHMARK_SAMPLES) | grep AVERAGE | tee -a benchmark.out

benchmark-dsyevd: full
> ./zig-out/x86_64-linux-musl/benchmark/dsyevd 2000 $(BENCHMARK_SAMPLES) | grep AVERAGE | tee -a benchmark.out

# DOCUMENTATION TARGETS ================================================================================================================================================================================

docs: library
> zig build docs && cp -r education/python docs && ./script/docstopdf.sh && cd docs && bundler install

serve: docs
> cd docs && bundle exec jekyll serve

# ANALYSIS TARGETS =====================================================================================================================================================================================

lines:
> git ls-files | grep -E '\.cpp$$|\.h$$|\.zig$$' | xargs cat | awk 'NF' | wc -l

linguist:
> github-linguist

# LIBRARY TARGETS ======================================================================================================================================================================================

library: external-x86_64-linux/.done

external-x86_64-linux:
> mkdir -p external-x86_64-linux/include external-x86_64-linux/lib

external-x86_64-linux/.done: zig-bin/.done external-x86_64-linux/eigen.tar.gz external-x86_64-linux/libint.tar.gz external-x86_64-linux/llvm.tar.xz external-x86_64-linux/openblas.tar.gz external-x86_64-linux/include/exprtk.hpp
> ./script/library.sh && touch $@

external-x86_64-linux/eigen.tar.gz: | external-x86_64-linux
> wget -q -O $@ https://gitlab.com/libeigen/eigen/-/archive/$(EIGEN_VERSION)/eigen-$(EIGEN_VERSION).tar.gz

external-x86_64-linux/libint.tar.gz: | external-x86_64-linux
> wget -q -O $@ https://github.com/evaleev/libint/releases/download/v$(LIBINT_VERSION)/libint-$(LIBINT_VERSION).tgz

external-x86_64-linux/llvm.tar.xz: | external-x86_64-linux
> wget -q -O $@ https://github.com/llvm/llvm-project/releases/download/llvmorg-$(LLVM_VERSION)/llvm-project-$(LLVM_VERSION).src.tar.xz

external-x86_64-linux/openblas.tar.gz: | external-x86_64-linux
> wget -q -O $@ https://github.com/OpenMathLib/OpenBLAS/releases/download/v$(OPENBLAS_VERSION)/OpenBLAS-$(OPENBLAS_VERSION).tar.gz

external-x86_64-linux/include/exprtk.hpp: | external-x86_64-linux
> wget -q -O $@ https://raw.githubusercontent.com/ArashPartow/exprtk/refs/heads/master/exprtk.hpp


# COMPILER TARGETS =====================================================================================================================================================================================

zig-bin:
> mkdir -p zig-bin

zig-bin/.done: zig-bin/zig zig-bin/zls
> touch $@

zig-bin/zig: | zig-bin
> wget -q -O zig.tar.xz https://ziglang.org/download/$(ZIG_VERSION)/zig-x86_64-linux-$(ZIG_VERSION).tar.xz && tar -xf zig.tar.xz -C zig-bin --strip-components=1 && rm zig.tar.xz

zig-bin/zls: | zig-bin
> wget -q -O zls.tar.xz https://github.com/zigtools/zls/releases/download/$(ZLS_VERSION)/zls-x86_64-linux.tar.xz && tar -xf zls.tar.xz -C zig-bin && rm zls.tar.xz

# CLEAN TARGETS ========================================================================================================================================================================================

clean:
> git clean -dffx

clean-cache:
> rm -rf zig-cache zig-global

clean-output:
> rm -rf zig-out
