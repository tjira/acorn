const std = @import("std"); const builtin = @import("builtin");

const targets: []const std.Target.Query = &.{
    .{.os_tag = .linux  , .cpu_arch = .aarch64},
    .{.os_tag = .linux  , .cpu_arch = .arm    },
    .{.os_tag = .linux  , .cpu_arch = .riscv64},
    .{.os_tag = .linux  , .cpu_arch = .x86    },
    .{.os_tag = .linux  , .cpu_arch = .x86_64 },
    .{.os_tag = .macos  , .cpu_arch = .aarch64},
    .{.os_tag = .macos  , .cpu_arch = .x86_64 },
    .{.os_tag = .wasi   , .cpu_arch = .wasm32 },
    .{.os_tag = .wasi   , .cpu_arch = .wasm64 },
    .{.os_tag = .windows, .cpu_arch = .aarch64},
    .{.os_tag = .windows, .cpu_arch = .x86    },
    .{.os_tag = .windows, .cpu_arch = .x86_64 },
};

pub fn build(builder: *std.Build) !void {
    const debug = builder.option(bool, "DEBUG", "Build everything in the debug mode") orelse false;

    for (targets) |target| {

        const main_executable = builder.addExecutable(.{
            .name = "acorn",
            .optimize = if (debug) .Debug else .ReleaseFast,
            .root_source_file = builder.path("src/main.zig"),
            .single_threaded = true,
            .strip = !debug,
            .target = builder.resolveTargetQuery(target)
        });

        const main_install = builder.addInstallArtifact(main_executable, .{
            .dest_dir = .{.override = .{.custom = try target.zigTriple(builder.allocator)}}
        });

        builder.getInstallStep().dependOn(&main_install.step);

        if (builtin.target.cpu.arch == target.cpu_arch and builtin.target.os.tag == target.os_tag) {

            const docs_install = builder.addInstallDirectory(.{
                .install_dir = .prefix, .install_subdir = "docs", .source_dir = builder.addExecutable(.{
                    .name = "main", .root_source_file = builder.path("src/main.zig"), .target = builder.host
                }).getEmittedDocs(),
            });

            var script_step = builder.step("script", "Generate executable scripts"); script_step.makeFn = script; script_step.dependOn(builder.getInstallStep());

            builder.step("docs", "Compile code documentation" ).dependOn(&docs_install                           .step);
            builder.step("run",  "Run the compiled executable").dependOn(&builder.addRunArtifact(main_executable).step);
        }
    }

    const test_executable = builder.addTest(.{
        .name = "test",
        .optimize = if (debug) .Debug else .ReleaseFast,
        .root_source_file = builder.path("test/main.zig"),
        .single_threaded = true,
        .strip = !debug,
        .target = builder.host
    });

    test_executable.root_module.addImport("acorn", builder.addModule("acorn", .{.root_source_file = builder.path("src/main.zig")}));

    builder.step("test", "Run unit tests").dependOn(&builder.addRunArtifact(test_executable).step);
}

pub fn script(self: *std.Build.Step, progress: std.Progress.Node) !void {
    _ = self; _ = progress;

    for (targets) |target| if (target.os_tag == .linux) {
        try linuxScripts(std.heap.page_allocator, try target.zigTriple(std.heap.page_allocator));
    };
}

pub fn linuxScripts(allocator: std.mem.Allocator, target: []const u8) !void {
    const path = try std.mem.concat(allocator, u8, &[_][]const u8{"zig-out/", target, "/"});

    const mode = switch (builtin.os.tag) {
        .wasi => 0, .windows => 0, else => 0o755
    };

    const file_matsort  = try std.fs.cwd().createFile(try std.mem.concat(allocator, u8, &[_][]const u8{path, "matsort" }), .{.mode=mode});
    const file_mersenne = try std.fs.cwd().createFile(try std.mem.concat(allocator, u8, &[_][]const u8{path, "mersenne"}), .{.mode=mode});
    const file_prime    = try std.fs.cwd().createFile(try std.mem.concat(allocator, u8, &[_][]const u8{path, "prime"   }), .{.mode=mode});

    const header =
        \\#!/usr/bin/bash
        \\
        \\clean() {{ rm -f input.json; }}; trap clean SIGINT
    ;

    const content_matsort =
        \\usage() {{ cat <<EOF | sed 's/^[[:space:]]*//'
        \\  Usage: $(basename "$0") [options]
        \\
        \\  Options:
        \\    -m <count> Matrix to sort.
        \\    -h         Display this help message and exit.
        \\EOF
        \\}}
        \\
        \\while getopts "m:h" OPT; do
        \\  case "$OPT" in
        \\    m ) MATRIX="$OPTARG" ;;
        \\    h ) usage && exit 0 ;;
        \\  esac
        \\done
        \\
        \\echo '{{
        \\  "sort" : {{
        \\    "input" : "'"$MATRIX"'",
        \\    "algorithm" : "bubble",
        \\    "output" : "'"$MATRIX"'"
        \\  }}
        \\}}' > input.json
    ;

    const content_mersenne =
        \\COUNT=10; OUTPUT="mersenne.mat"; START=1
        \\
        \\usage() {{ cat <<EOF | sed 's/^[[:space:]]*//'
        \\  Usage: $(basename "$0") [options]
        \\
        \\  Options:
        \\    -c <count>  Number of Mersenne primes to generate. (default: ${{COUNT}})
        \\    -o <output> Output file. (default: ${{OUTPUT}})
        \\    -s <start>  Starting number. (default: ${{START}})
        \\    -h          Display this help message and exit.
        \\EOF
        \\}}
        \\
        \\while getopts "c:o:s:h" OPT; do
        \\  case "$OPT" in
        \\    c ) COUNT="$OPTARG" ;;
        \\    o ) OUTPUT="$OPTARG" ;;
        \\    s ) START="$OPTARG" ;;
        \\    h ) usage && exit 0 ;;
        \\  esac
        \\done
        \\
        \\echo '{{
        \\  "prime" : {{
        \\    "mode" : "mersenne",
        \\    "generate" : {{
        \\      "count" : '"$COUNT"',
        \\      "output" : "'"$OUTPUT"'"
        \\    }},
        \\    "number" : '"$START"'
        \\  }}
        \\}}' > input.json
    ;

    const content_prime =
        \\COUNT=10; LOG_INTERVAL=1; OUTPUT="primes.mat"; START=1
        \\
        \\usage() {{ cat <<EOF | sed 's/^[[:space:]]*//'
        \\   Usage: $(basename "$0") [options]
        \\
        \\   Options:
        \\     -c <count>        Number of primes to generate. (default: ${{COUNT}})
        \\     -l <log_interval> Log interval for output. (default: ${{LOG_INTERVAL}})
        \\     -o <output>       Output file. (default: ${{OUTPUT}})
        \\     -s <start>        Starting number. (default: ${{START}})
        \\     -h                Display this help message and exit.
        \\EOF
        \\}}
        \\
        \\while getopts "c:l:o:s:h" OPT; do
        \\  case "$OPT" in
        \\    c ) COUNT="$OPTARG" ;;
        \\    l ) LOG_INTERVAL="$OPTARG" ;;
        \\    o ) OUTPUT="$OPTARG" ;;
        \\    s ) START="$OPTARG" ;;
        \\    h ) usage && exit 0 ;;
        \\  esac
        \\done
        \\
        \\echo '{{
        \\  "prime" : {{
        \\    "mode" : "basic",
        \\    "generate" : {{
        \\      "count" : '"$COUNT"',
        \\      "log_interval" : '"$LOG_INTERVAL"',
        \\      "output" : "'"$OUTPUT"'"
        \\    }},
        \\    "number" : '"$START"'
        \\  }}
        \\}}' > input.json
    ;

    try file_matsort .writer().print(header ++ "\n\n" ++ content_matsort  ++ "\n\nacorn && clean", .{}); file_matsort .close();
    try file_mersenne.writer().print(header ++ "\n\n" ++ content_mersenne ++ "\n\nacorn && clean", .{}); file_mersenne.close();
    try file_prime   .writer().print(header ++ "\n\n" ++ content_prime    ++ "\n\nacorn && clean", .{});    file_prime.close();
}
