.RECIPEPREFIX = >

ZIG_VERSION = 0.14.1
ZLS_VERSION = 0.14.0

all: acorn benchmark

# ACORN TARGETS ========================================================================================================================================================================================

acorn: library
> ./zig-bin/zig build

benchmark: library
> ./zig-bin/zig build -DBENCHMARK

run: acorn
> ./zig-bin/zig build run

test: library
> ./zig-bin/zig build test

# DOCUMENTATION TARGETS ================================================================================================================================================================================

serve: compiler
> ./zig-bin/zig build-obj -Iinclude -Iexternal-x86_64-linux/include -lc -femit-docs=docs/code src/main.zig && rm main.o* && python -m http.server 8000 -d docs/code

# LIBRARY TARGETS ======================================================================================================================================================================================

library: compiler
> [ ! -d external-x86_64-linux ] && ./script/library.sh || true

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
