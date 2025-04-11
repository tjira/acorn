const std = @import("std"); const builtin = @import("builtin");

const targets: []const std.Target.Query = &.{
    .{.os_tag = .linux  , .cpu_arch = .aarch64},
    .{.os_tag = .linux  , .cpu_arch = .x86_64 },
    .{.os_tag = .macos  , .cpu_arch = .aarch64},
    .{.os_tag = .macos  , .cpu_arch = .x86_64 },
    .{.os_tag = .windows, .cpu_arch = .aarch64},
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

        main_executable.addIncludePath(.{.cwd_relative="/usr/include"});
        main_executable.addLibraryPath(.{.cwd_relative="/usr/lib"    });

        main_executable.linkLibC();
        // main_executable.linkSystemLibrary("lapacke");

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

    for (targets) |target| if (target.os_tag == .windows) {
        try windowsScripts(std.heap.page_allocator, try target.zigTriple(std.heap.page_allocator));
    };
}

pub fn linuxScripts(allocator: std.mem.Allocator, target: []const u8) !void {
    const path = try std.mem.concat(allocator, u8, &[_][]const u8{"zig-out/", target, "/"});

    const mode = switch (builtin.os.tag) {
        .wasi => 0, .windows => 0, else => 0o755
    };

    const file_fibonacci = try std.fs.cwd().createFile(try std.mem.concat(allocator, u8, &[_][]const u8{path, "fibonacci"}), .{.mode=mode});
    const file_matsort   = try std.fs.cwd().createFile(try std.mem.concat(allocator, u8, &[_][]const u8{path, "matsort"  }), .{.mode=mode});
    const file_mersenne  = try std.fs.cwd().createFile(try std.mem.concat(allocator, u8, &[_][]const u8{path, "mersenne" }), .{.mode=mode});
    const file_prime     = try std.fs.cwd().createFile(try std.mem.concat(allocator, u8, &[_][]const u8{path, "prime"    }), .{.mode=mode});

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

pub fn windowsScripts(allocator: std.mem.Allocator, target: []const u8) !void {
    const path = try std.mem.concat(allocator, u8, &[_][]const u8{"zig-out/", target, "/"});

    const mode = switch (builtin.os.tag) {
        .wasi => 0, .windows => 0, else => 0o755
    };

    const file_matsort  = try std.fs.cwd().createFile(try std.mem.concat(allocator, u8, &[_][]const u8{path, "matsort.ps1" }), .{.mode=mode});
    const file_mersenne = try std.fs.cwd().createFile(try std.mem.concat(allocator, u8, &[_][]const u8{path, "mersenne.ps1"}), .{.mode=mode});
    const file_prime    = try std.fs.cwd().createFile(try std.mem.concat(allocator, u8, &[_][]const u8{path, "prime.ps1"   }), .{.mode=mode});

    const header =
        \\$IDX = 0
        \\
        \\function Show-Usage {{
        \\  $NAME = [System.IO.Path]::GetFileNameWithoutExtension($PSCommandPath)
        \\  $USAGE = @"
        \\Usage: $NAME [options]
        \\
        \\Options:
    ;

    const footer =
        \\[System.IO.File]::WriteAllText(".input.json", $CONTENT, [System.Text.UTF8Encoding]::new($false))
        \\
        \\& acorn .input.json
        \\
        \\if ($LASTEXITCODE -eq 0) {{
        \\  Remove-Item -Force input.json -ErrorAction SilentlyContinue
        \\}}
    ;

    const content_matsort =
        \\  -m <start> Matrix to sort.
        \\  -h         Display this help message and exit.
        \\"@
        \\  Write-Output $USAGE.Trim()
        \\}}
        \\
        \\while ($IDX -lt $args.Count) {{ switch ($args[$IDX]) {{
        \\  '-m' {{ $IDX++; if ($IDX -lt $args.Count) {{ $MATRIX = $args[$IDX] }} }}
        \\  '-h' {{ Show-Usage; exit 0 }} Default {{ Show-Usage; exit 1 }}
        \\}} $IDX++ }}
        \\
        \\if ($MATRIX -ne "null") {{ $MATRIX = '"' + $MATRIX + '"' }}
        \\
        \\$CONTENT = "{{`n" +
        \\"  `"sort`" : {{`n" +
        \\"    `"input`" : $MATRIX,`n" +
        \\"    `"algorithm`" : `"bubble`",`n" +
        \\"    `"output`" : $MATRIX`n" +
        \\"  }}`n" +
        \\"}}"
    ;

    const content_mersenne =
        \\  -c <count>  Number of Mersenne primes to generate. (default: $COUNT)
        \\  -o <output> Output file. (default: $OUTPUT)
        \\  -s <start>  Starting number. (default: $START)
        \\  -h          Display this help message and exit.
        \\"@
        \\  Write-Output $USAGE.Trim()
        \\}}
        \\
        \\$COUNT = 10; $OUTPUT = "null"; $START = 1
        \\
        \\while ($IDX -lt $args.Count) {{ switch ($args[$IDX]) {{
        \\  '-c' {{ $IDX++; if ($IDX -lt $args.Count) {{ $COUNT = $args[$IDX] }} }}
        \\  '-o' {{ $IDX++; if ($IDX -lt $args.Count) {{ $OUTPUT = $args[$IDX] }} }}
        \\  '-s' {{ $IDX++; if ($IDX -lt $args.Count) {{ $START = $args[$IDX] }} }}
        \\  '-h' {{ Show-Usage; exit 0 }} Default {{ Show-Usage; exit 1 }}
        \\}} $IDX++ }}
        \\
        \\if ($OUTPUT -ne "null") {{ $OUTPUT = '"' + $OUTPUT + '"' }}
        \\
        \\$CONTENT = "{{`n" +
        \\"  `"prime`" : {{`n" +
        \\"    `"mode`" : `"mersenne`",`n" +
        \\"    `"generate`" : {{`n" +
        \\"      `"count`" : $COUNT,`n" +
        \\"      `"output`" : $OUTPUT`n" +
        \\"    }},`n" +
        \\"    `"number`" : $START`n" +
        \\"  }}`n" +
        \\"}}"
    ;

    const content_prime =
        \\  -c <count>        Number of primes to generate. (default: $COUNT)
        \\  -l <log_interval> Log interval for output. (default: $LOG_INTERVAL)
        \\  -o <output>       Output file. (default: $OUTPUT)
        \\  -s <start>        Starting number. (default: $START)
        \\  -h                Display this help message and exit.
        \\"@
        \\  Write-Output $USAGE.Trim()
        \\}}
        \\
        \\$COUNT = 10; $LOG_INTERVAL = 1; $OUTPUT = "null"; $START = 1
        \\
        \\while ($IDX -lt $args.Count) {{ switch ($args[$IDX]) {{
        \\  '-c' {{ $IDX++; if ($IDX -lt $args.Count) {{ $COUNT = $args[$IDX] }} }}
        \\  '-l' {{ $IDX++; if ($IDX -lt $args.Count) {{ $LOG_INTERVAL = $args[$IDX] }} }}
        \\  '-o' {{ $IDX++; if ($IDX -lt $args.Count) {{ $OUTPUT = $args[$IDX] }} }}
        \\  '-s' {{ $IDX++; if ($IDX -lt $args.Count) {{ $START = $args[$IDX] }} }}
        \\  '-h' {{ Show-Usage; exit 0 }} Default {{ Show-Usage; exit 1 }}
        \\}} $IDX++ }}
        \\
        \\if ($OUTPUT -ne "null") {{ $OUTPUT = '"' + $OUTPUT + '"' }}
        \\
        \\$CONTENT = "{{`n" +
        \\"  `"prime`" : {{`n" +
        \\"    `"mode`" : `"basic`",`n" +
        \\"    `"generate`" : {{`n" +
        \\"      `"count`" : $COUNT,`n" +
        \\"      `"log_interval`" : $LOG_INTERVAL,`n" +
        \\"      `"output`" : $OUTPUT`n" +
        \\"    }},`n" +
        \\"    `"number`" : $START`n" +
        \\"  }}`n" +
        \\"}}"
    ;

    try file_matsort .writer().print(header ++ "\n\n" ++ content_matsort  ++ "\n\n" ++ footer, .{}); file_matsort .close();
    try file_mersenne.writer().print(header ++ "\n\n" ++ content_mersenne ++ "\n\n" ++ footer, .{}); file_mersenne.close();
    try file_prime   .writer().print(header ++ "\n\n" ++ content_prime    ++ "\n\n" ++ footer, .{});    file_prime.close();
}
