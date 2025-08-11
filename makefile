all: library acorn script

.PHONY: acorn run script test

acorn:
	./zig-bin/zig build

library:
	[ ! -d external-x86_64-linux ] && ./script/library.sh || true

run:
	./zig-bin/zig build run

script:
	./zig-bin/zig build script

test:
	./zig-bin/zig build test
