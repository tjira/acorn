const std = @import("std"); const builtin = @import("builtin");

const target: std.Target.Query = .{.os_tag = builtin.target.os.tag, .cpu_arch = builtin.target.cpu.arch, .abi = .gnu};

pub fn build(builder: *std.Build) !void {
    const debug = builder.option(bool, "DEBUG", "Build everything in the debug mode") orelse false;

    const main_executable = builder.addExecutable(.{
        .name = "acorn",
        .optimize = if (debug) .Debug else .ReleaseFast,
        .root_source_file = builder.path("src/main.zig"),
        .strip = !debug,
        .target = builder.resolveTargetQuery(target)
    });

    const test_executable = builder.addTest(.{
        .name = "test",
        .optimize = if (debug) .Debug else .ReleaseFast,
        .root_source_file = builder.path("test/main.zig"),
        .strip = !debug,
        .target = builder.host
    });

    main_executable.addIncludePath(.{.cwd_relative="/usr/include"}); test_executable.addIncludePath(.{.cwd_relative="/usr/include"});
    main_executable.addLibraryPath(.{.cwd_relative="/usr/lib"    }); test_executable.addLibraryPath(.{.cwd_relative="/usr/lib"    });

    main_executable.linkLibC(); test_executable.linkLibC();

    main_executable.linkSystemLibrary("cblas"  ); test_executable.linkSystemLibrary("cblas"  );
    main_executable.linkSystemLibrary("fftw3"  ); test_executable.linkSystemLibrary("fftw3"  );
    main_executable.linkSystemLibrary("gsl"    ); test_executable.linkSystemLibrary("gsl"    );
    main_executable.linkSystemLibrary("lapacke"); test_executable.linkSystemLibrary("lapacke");

    test_executable.root_module.addImport("acorn", builder.addModule("acorn", .{.root_source_file = builder.path("src/main.zig")}));

    const docs_install = builder.addInstallDirectory(.{
        .install_dir = .prefix, .install_subdir = "docs", .source_dir = main_executable.getEmittedDocs()
    });

    const main_install = builder.addInstallArtifact(main_executable, .{
        .dest_dir = .{.override = .{.custom = try target.zigTriple(builder.allocator)}}
    });

    builder.getInstallStep().dependOn(&main_install.step);

    builder.step("docs", "Compile code documentation" ).dependOn(&docs_install                           .step);
    builder.step("run",  "Run the compiled executable").dependOn(&builder.addRunArtifact(main_executable).step);
    builder.step("test", "Run unit tests"             ).dependOn(&builder.addRunArtifact(test_executable).step);

    var script_step = builder.step("script", "Generate executable scripts"); script_step.makeFn = script; script_step.dependOn(builder.getInstallStep());
}

pub fn script(self: *std.Build.Step, progress: std.Progress.Node) !void {
    _ = self; _ = progress;

    if (target.os_tag == .linux) {
        try linuxScripts(std.heap.page_allocator, try target.zigTriple(std.heap.page_allocator));
    }
}

pub fn linuxScripts(allocator: std.mem.Allocator, bindir: []const u8) !void {
    const path = try std.mem.concat(allocator, u8, &[_][]const u8{"zig-out/", bindir, "/"});

    const file_fibonacci = try std.fs.cwd().createFile(try std.mem.concat(allocator, u8, &[_][]const u8{path, "fibonacci"}), .{.mode=0o755});
    const file_matsort   = try std.fs.cwd().createFile(try std.mem.concat(allocator, u8, &[_][]const u8{path, "matsort"  }), .{.mode=0o755});
    const file_mersenne  = try std.fs.cwd().createFile(try std.mem.concat(allocator, u8, &[_][]const u8{path, "mersenne" }), .{.mode=0o755});
    const file_prime     = try std.fs.cwd().createFile(try std.mem.concat(allocator, u8, &[_][]const u8{path, "prime"    }), .{.mode=0o755});

    const header =
        \\#!/usr/bin/bash
        \\
        \\clean() {{ rm -f .input.json; }}; trap clean SIGINT
        \\
        \\usage() {{
        \\  cat <<EOF
        \\Usage: $(basename $0) [options]
        \\
        \\Options:
    ;

    const content_fibonacci =
        \\  -c <count>        Number of Fibonacci numbers to generate. (default: ${{COUNT}})
        \\  -l <log_interval> Log interval for output. (default: ${{LOG_INTERVAL}})
        \\  -o <output>       Output file. (default: ${{OUTPUT}})
        \\  -h                Display this help message and exit.
        \\EOF
        \\
        \\}}; COUNT=10; LOG_INTERVAL=1; OUTPUT=null
        \\
        \\while getopts "c:l:o:h" OPT; do case "$OPT" in
        \\  c ) COUNT="$OPTARG" ;; l ) LOG_INTERVAL="$OPTARG" ;; o ) OUTPUT="$OPTARG" ;; h ) usage && exit 0 ;; \? ) usage && exit 1 ;;
        \\esac done
        \\
        \\[[ "$OUTPUT" != "null" ]] && OUTPUT="\"$OUTPUT\""
        \\
        \\echo '{{
        \\  "fibonacci" : {{
        \\    "count" : '"$COUNT"',
        \\    "log_interval" : '"$LOG_INTERVAL"',
        \\    "output" : '"$OUTPUT"'
        \\  }}
        \\}}' > .input.json
    ;

    const content_matsort =
        \\  -m <count> Matrix to sort.
        \\  -h         Display this help message and exit.
        \\EOF
        \\
        \\}}
        \\
        \\while getopts "m:h" OPT; do case "$OPT" in
        \\  m ) MATRIX="$OPTARG" ;; h ) usage && exit 0 ;; \? ) usage && exit 1 ;;
        \\esac done
        \\
        \\echo '{{
        \\  "sort" : {{
        \\    "input" : "'"$MATRIX"'",
        \\    "algorithm" : "bubble",
        \\    "output" : "'"$MATRIX"'"
        \\  }}
        \\}}' > .input.json
    ;

    const content_mersenne =
        \\  -c <count>  Number of Mersenne primes to generate. (default: ${{COUNT}})
        \\  -o <output> Output file. (default: ${{OUTPUT}})
        \\  -s <start>  Starting number. (default: ${{START}})
        \\  -h          Display this help message and exit.
        \\EOF
        \\
        \\}}; COUNT=10; OUTPUT=null; START=1
        \\
        \\while getopts "c:o:s:h" OPT; do case "$OPT" in
        \\  c ) COUNT="$OPTARG" ;; o ) OUTPUT="$OPTARG" ;; s ) START="$OPTARG" ;; h ) usage && exit 0 ;; \? ) usage && exit 1 ;;
        \\esac done
        \\
        \\[[ "$OUTPUT" != "null" ]] && OUTPUT="\"$OUTPUT\""
        \\
        \\echo '{{
        \\  "prime" : {{
        \\    "mode" : "mersenne",
        \\    "generate" : {{
        \\      "count" : '"$COUNT"',
        \\      "output" : '"$OUTPUT"'
        \\    }},
        \\    "number" : '"$START"'
        \\  }}
        \\}}' > .input.json
    ;

    const content_prime =
        \\  -c <count>        Number of primes to generate. (default: ${{COUNT}})
        \\  -l <log_interval> Log interval for output. (default: ${{LOG_INTERVAL}})
        \\  -o <output>       Output file. (default: ${{OUTPUT}})
        \\  -s <start>        Starting number. (default: ${{START}})
        \\  -h                Display this help message and exit.
        \\EOF
        \\
        \\}}; COUNT=10; LOG_INTERVAL=1; OUTPUT=null; START=1
        \\
        \\while getopts "c:l:o:s:h" OPT; do case "$OPT" in
        \\  c ) COUNT="$OPTARG" ;; l ) LOG_INTERVAL="$OPTARG" ;; o ) OUTPUT="$OPTARG" ;; s ) START="$OPTARG" ;; h ) usage && exit 0 ;; \? ) usage && exit 1 ;;
        \\esac done
        \\
        \\[[ "$OUTPUT" != "null" ]] && OUTPUT="\"$OUTPUT\""
        \\
        \\echo '{{
        \\  "prime" : {{
        \\    "mode" : "basic",
        \\    "generate" : {{
        \\      "count" : '"$COUNT"',
        \\      "log_interval" : '"$LOG_INTERVAL"',
        \\      "output" : '"$OUTPUT"'
        \\    }},
        \\    "number" : '"$START"'
        \\  }}
        \\}}' > .input.json
    ;

    try file_fibonacci.writer().print(header ++ "\n" ++ content_fibonacci ++ "\n\nacorn .input.json && clean", .{}); file_fibonacci.close();
    try file_matsort  .writer().print(header ++ "\n" ++ content_matsort   ++ "\n\nacorn .input.json && clean", .{});  file_matsort .close();
    try file_mersenne .writer().print(header ++ "\n" ++ content_mersenne  ++ "\n\nacorn .input.json && clean", .{});  file_mersenne.close();
    try file_prime    .writer().print(header ++ "\n" ++ content_prime     ++ "\n\nacorn .input.json && clean", .{});     file_prime.close();
}
