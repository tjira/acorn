.RECIPEPREFIX = >

ARCH := $(shell uname -m | tr '[:upper:]' '[:lower:]')
OS   := $(shell uname -s | tr '[:upper:]' '[:lower:]')

ZIG_VERSION := 0.15.1
ZLS_VERSION := 0.15.0

all: acorn

# ACORN BUILDING TARGETS ===============================================================================================================================================================================

.PHONY: acorn

acorn: zig-bin/zig
> ./zig-bin/zig build

# EXTERNAL TARGETS =====================================================================================================================================================================================

zig-bin/zig:
> mkdir -p zig-bin && wget -q -O - https://ziglang.org/download/$(ZIG_VERSION)/zig-$(ARCH)-$(OS)-$(ZIG_VERSION).tar.xz       | tar -Jx -C zig-bin --strip-components=1 && touch zig-bin/zig
> mkdir -p zig-bin && wget -q -O - https://github.com/zigtools/zls/releases/download/$(ZLS_VERSION)/zls-$(ARCH)-$(OS).tar.xz | tar -Jx -C zig-bin --strip-components=0 && touch zig-bin/zls

# CLEAN TARGETS ========================================================================================================================================================================================

clean:
> git clean -dffx
