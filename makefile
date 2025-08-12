.RECIPEPREFIX = >

ZIG_VERSION = 0.14.1
ZLS_VERSION = 0.14.0

all: acorn benchmark

# ACORN TARGETS ========================================================================================================================================================================================

acorn: library
> ./zig-bin/zig build

test: library
> ./zig-bin/zig build test

# DOCUMENTATION TARGETS ================================================================================================================================================================================

docs: library
> ./zig-bin/zig build-obj -Iinclude -Iexternal-x86_64-linux/include -lc -femit-docs=code src/main.zig && rm *.o && cp -r education/python docs && ./script/docstopdf.sh && cd docs && bundler install

serve: docs
> cd docs && bundle exec jekyll serve

# LIBRARY TARGETS ======================================================================================================================================================================================

library: compiler
> [ ! -d external-x86_64-linux ] && wget -q -O boost.tar.gz    https://github.com/boostorg/boost/releases/download/boost-1.88.0/boost-1.88.0-b2-nodocs.tar.gz || true
> [ ! -d external-x86_64-linux ] && wget -q -O eigen.tar.gz                              https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz || true
> [ ! -d external-x86_64-linux ] && wget -q -O libint.tar.gz                    https://github.com/evaleev/libint/releases/download/v2.11.1/libint-2.11.1.tgz || true
> [ ! -d external-x86_64-linux ] && wget -q -O llvm.tar.gz                       https://github.com/llvm/llvm-project/archive/refs/tags/llvmorg-20.1.8.tar.gz || true
> [ ! -d external-x86_64-linux ] && wget -q -O openblas.tar.gz       https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.30/OpenBLAS-0.3.30.tar.gz || true
> [ ! -d external-x86_64-linux ] && ./script/library.sh && rm *.tar.gz || true

# COMPILER TARGETS =====================================================================================================================================================================================

compiler: zig-bin/zig zig-bin/zls

zig-bin:
> mkdir -p zig-bin

zig-bin/zig: | zig-bin
> wget -q -O zig.tar.xz https://ziglang.org/download/$(ZIG_VERSION)/zig-x86_64-linux-$(ZIG_VERSION).tar.xz
> tar -xf zig.tar.xz -C zig-bin --strip-components=1
> rm zig.tar.xz

zig-bin/zls: | zig-bin
> wget -q -O zls.tar.xz https://github.com/zigtools/zls/releases/download/$(ZLS_VERSION)/zls-x86_64-linux.tar.xz
> tar -xf zls.tar.xz -C zig-bin
> rm zls.tar.xz
