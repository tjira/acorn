const std = @import("std"); const builtin = @import("builtin");

const targets: []const std.Target.Query = &.{
    .{.os_tag = .linux, .cpu_arch = .x86_64, .abi = .gnu },
    .{.os_tag = .linux, .cpu_arch = .x86_64, .abi = .musl},
};

pub fn build(builder: *std.Build) !void {
    const debug = builder.option(bool, "DEBUG", "Build everything in the debug mode") orelse false;

    for (targets) |target| {

        const os_name   = switch (target.os_tag.?  ) {.linux  => "linux",  .macos   => "macos",   .windows => "windows", else => unreachable};
        const arch_name = switch (target.cpu_arch.?) {.x86_64 => "x86_64", .aarch64 => "aarch64", .riscv64 => "riscv64", else => unreachable};

        const benchmark_executable = builder.addExecutable(.{
            .name = "benchmark",
            .root_module = builder.createModule(.{
                .optimize = if (debug) .Debug else .ReleaseFast,
                .root_source_file = builder.path("benchmark/main.zig"),
                .strip = !debug,
                .target = builder.resolveTargetQuery(target)
            }),
            .use_llvm = true
        });

        const main_executable = builder.addExecutable(.{
            .name = "acorn",
            .root_module = builder.createModule(.{
                .optimize = if (debug) .Debug else .ReleaseFast,
                .root_source_file = builder.path("src/main.zig"),
                .strip = !debug,
                .target = builder.resolveTargetQuery(target)
            }),
            .use_llvm = true
        });

        const test_executable = builder.addTest(.{
            .name = "test",
            .root_module = builder.createModule(.{
                .optimize = if (debug) .Debug else .ReleaseFast,
                .root_source_file = builder.path("test/main.zig"),
                .strip = !debug,
                .target = builder.resolveTargetQuery(target)
            }),
            .use_llvm = true
        });

        const external_include = try std.mem.concat(builder.allocator, u8, &[_][]const u8{"external-", arch_name, "-", os_name, "/include"});
        const external_lib     = try std.mem.concat(builder.allocator, u8, &[_][]const u8{"external-", arch_name, "-", os_name, "/lib"    });

        main_executable.root_module.addIncludePath(.{.cwd_relative = "include"       });
        main_executable.root_module.addIncludePath(.{.cwd_relative = external_include});
        main_executable.root_module.addLibraryPath(.{.cwd_relative = external_lib    });

        main_executable.root_module.addIncludePath(.{.cwd_relative = try std.mem.concat(builder.allocator, u8, &[_][]const u8{external_include, "/eigen3"})});

        main_executable.root_module.addCSourceFile(.{.file = builder.path("src/eigen.cpp" ), .flags = &[_][]const u8{"-fopenmp"}, .language = .cpp});
        main_executable.root_module.addCSourceFile(.{.file = builder.path("src/exprtk.cpp"), .flags = &[_][]const u8{"-fopenmp"}, .language = .cpp});
        main_executable.root_module.addCSourceFile(.{.file = builder.path("src/libint.cpp"), .flags = &[_][]const u8{"-fopenmp"}, .language = .cpp});

        main_executable.root_module.link_libc   = true;
        main_executable.root_module.link_libcpp = true;

        main_executable.root_module.linkSystemLibrary("int2",     .{.preferred_link_mode = .static});
        main_executable.root_module.linkSystemLibrary("omp",      .{.preferred_link_mode = .static});
        main_executable.root_module.linkSystemLibrary("openblas", .{.preferred_link_mode = .static});

        test_executable     .root_module.addImport("acorn", main_executable.root_module);
        benchmark_executable.root_module.addImport("acorn", main_executable.root_module);

        const main_executable_install = builder.addInstallArtifact(main_executable, .{
            .dest_dir = .{.override = .{.custom = try target.zigTriple(builder.allocator)}}
        });

        const benchmark_executable_install = builder.addInstallArtifact(benchmark_executable, .{
            .dest_dir = .{.override = .{.custom = try target.zigTriple(builder.allocator)}}
        });

        builder.getInstallStep().dependOn(     &main_executable_install.step);
        builder.getInstallStep().dependOn(&benchmark_executable_install.step);

        if (builtin.target.cpu.arch == target.cpu_arch and builtin.target.os.tag == target.os_tag and target.abi == .musl) {

            builder.step("run",  "Run the compiled executable"   ).dependOn(&builder.addRunArtifact(main_executable).step);
            builder.step("test", "Run unit tests"                ).dependOn(&builder.addRunArtifact(test_executable).step);

            const docs = builder.addInstallDirectory(.{.source_dir = main_executable.getEmittedDocs(), .install_dir = .{.custom = "../docs"}, .install_subdir = "code"});

            builder.step("docs", "Install the code documentation").dependOn(&docs.step);
        }
    }

    std.fs.cwd().makeDir("zig-out") catch {};

    for (targets) |target| if (target.os_tag == .linux) {
        try linuxScripts(try target.zigTriple(std.heap.page_allocator), std.heap.page_allocator);
    };
}

pub fn linuxScripts(bindir: []const u8, allocator: std.mem.Allocator) !void {
    const path = try std.mem.concat(allocator, u8, &[_][]const u8{"zig-out/", bindir, "/"});

    std.fs.cwd().makeDir(path) catch {};

    const file_fibonacci = try std.fs.cwd().createFile(try std.mem.concat(allocator, u8, &[_][]const u8{path, "fibonacci"}), .{.mode=0o755});
    const file_hf        = try std.fs.cwd().createFile(try std.mem.concat(allocator, u8, &[_][]const u8{path, "hf"       }), .{.mode=0o755});
    const file_integral  = try std.fs.cwd().createFile(try std.mem.concat(allocator, u8, &[_][]const u8{path, "integral" }), .{.mode=0o755});
    const file_matmul    = try std.fs.cwd().createFile(try std.mem.concat(allocator, u8, &[_][]const u8{path, "matmul"   }), .{.mode=0o755});
    const file_matsort   = try std.fs.cwd().createFile(try std.mem.concat(allocator, u8, &[_][]const u8{path, "matsort"  }), .{.mode=0o755});
    const file_mersenne  = try std.fs.cwd().createFile(try std.mem.concat(allocator, u8, &[_][]const u8{path, "mersenne" }), .{.mode=0o755});
    const file_prime     = try std.fs.cwd().createFile(try std.mem.concat(allocator, u8, &[_][]const u8{path, "prime"    }), .{.mode=0o755});
    const file_randmat   = try std.fs.cwd().createFile(try std.mem.concat(allocator, u8, &[_][]const u8{path, "randmat"  }), .{.mode=0o755});
    const file_schlitter = try std.fs.cwd().createFile(try std.mem.concat(allocator, u8, &[_][]const u8{path, "schlitter"}), .{.mode=0o755});

    const header =
        \\#!/bin/bash
        \\
        \\INPUT_JSON=$(mktemp)
        \\
        \\usage() {
        \\  cat <<EOF
        \\Usage: $(basename $0) [options]
        \\
        \\Options:
    ;

    const content_fibonacci =
        \\  -c <count>        Number of Fibonacci numbers to generate. (default: ${COUNT})
        \\  -l <log_interval> Log interval for output. (default: ${LOG_INTERVAL})
        \\  -o <output>       Output file. (default: ${OUTPUT})
        \\  -h                Display this help message and exit.
        \\EOF
        \\
        \\}; COUNT=10; LOG_INTERVAL=1; OUTPUT=null
        \\
        \\while getopts "c:l:o:h" OPT; do case "$OPT" in
        \\  c ) COUNT="$OPTARG" ;;
        \\  l ) LOG_INTERVAL="$OPTARG" ;;
        \\  o ) OUTPUT="$OPTARG" ;;
        \\  h ) usage && exit 0 ;;
        \\  \?) usage && exit 1 ;;
        \\esac done
        \\
        \\[[ "$OUTPUT" != "null" ]] && OUTPUT="\"$OUTPUT\""
        \\
        \\echo '{
        \\  "fibonacci" : {
        \\    "count" : '"$COUNT"',
        \\    "log_interval" : '"$LOG_INTERVAL"',
        \\    "output" : '"$OUTPUT"'
        \\  }
        \\}' > $INPUT_JSON
    ;

    const content_hf =
        \\  -b <basis>       Atomic orbital basis. (default: ${BASIS})
        \\  -c <charge>      Charge of the system. (default: ${CHARGE})
        \\  -d <direct>      Boolean direct flag. (default: ${DIRECT})
        \\  -e <export>      Boolean export flag. (default: ${EXPORT})
        \\  -g <generalized> Boolean generalized flag. (default: ${GENERALIZED})
        \\  -i <input>       Input system file in .xyz format. (default: ${SYSTEM})
        \\  -n <nthread>     Number of threads to use. (default: ${N})
        \\  -h               Display this help message and exit.
        \\EOF
        \\
        \\}; N=1; BASIS="sto-3g"; CHARGE=0; DIRECT=false; EXPORT=false; GENERALIZED=false; INPUT="molecule.xyz";
        \\
        \\while getopts "b:c:degi:n:h" OPT; do case "$OPT" in
        \\  b ) BASIS="$OPTARG" ;;
        \\  c ) CHARGE="$OPTARG" ;;
        \\  d ) DIRECT=true ;;
        \\  e ) EXPORT=true ;;
        \\  g ) GENERALIZED=true ;;
        \\  i ) INPUT="$OPTARG" ;;
        \\  n ) N="$OPTARG" ;;
        \\  h ) usage && exit 0 ;;
        \\  \?) usage && exit 1 ;;
        \\esac done
        \\
        \\[[ "$EXPORT" == true ]] && D_AO="\"S_AO.mat\"" || D_AO="null"
        \\[[ "$EXPORT" == true ]] && C_AO="\"T_AO.mat\"" || C_AO="null"
        \\[[ "$EXPORT" == true ]] && F_AO="\"V_AO.mat\"" || F_AO="null"
        \\
        \\echo '{
        \\  "hartree_fock" : {
        \\      "system_file" : "'"$INPUT"'",
        \\      "system" : {
        \\          "charge" : 0
        \\      },
        \\      "integral" : {
        \\          "basis" : "'"$BASIS"'"
        \\      },
        \\      "write" : {
        \\          "coefficient" : '"$C_AO"',
        \\          "density" : '"$D_AO"',
        \\          "fock" : '"$F_AO"'
        \\      },
        \\      "threshold" : 1e-8,
        \\      "maxiter" : 1000,
        \\      "dsize" : 5,
        \\      "generalized" : '$GENERALIZED',
        \\      "direct" : '$DIRECT'
        \\  }
        \\}' > $INPUT_JSON
    ;

    const content_integral =
        \\  -b <basis>   Atomic orbital basis. (default: ${BASIS})
        \\  -e <export>  Boolean export flag. (default: ${EXPORT})
        \\  -i <input>   Input system file in .xyz format. (default: ${SYSTEM})
        \\  -n <nthread> Number of threads to use. (default: ${N})
        \\  -h           Display this help message and exit.
        \\EOF
        \\
        \\}; N=1; BASIS="sto-3g"; EXPORT=false; INPUT="molecule.xyz";
        \\
        \\while getopts "b:ei:n:h" OPT; do case "$OPT" in
        \\  b ) BASIS="$OPTARG" ;;
        \\  e ) EXPORT=true ;;
        \\  i ) INPUT="$OPTARG" ;;
        \\  n ) N="$OPTARG" ;;
        \\  h ) usage && exit 0 ;;
        \\  \?) usage && exit 1 ;;
        \\esac done
        \\
        \\[[ "$EXPORT" == true ]] && S_AO="\"S_AO.mat\"" || S_AO="null"
        \\[[ "$EXPORT" == true ]] && T_AO="\"T_AO.mat\"" || T_AO="null"
        \\[[ "$EXPORT" == true ]] && V_AO="\"V_AO.mat\"" || V_AO="null"
        \\[[ "$EXPORT" == true ]] && J_AO="\"J_AO.mat\"" || J_AO="null"
        \\
        \\echo '{
        \\  "atomic_integral" : {
        \\      "system_file" : "'"$INPUT"'",
        \\      "basis" : "'"$BASIS"'",
        \\      "write" : {
        \\          "overlap" : '"$S_AO"',
        \\          "kinetic" : '"$T_AO"',
        \\          "nuclear" : '"$V_AO"',
        \\          "coulomb" : '"$J_AO"'
        \\      }
        \\  }
        \\}' > $INPUT_JSON
    ;

    const content_matmul =
        \\  -l <left>    Left multiplicant.
        \\  -o <output>  Output matrix. (default: ${OUTPUT})
        \\  -p <print>   Boolean print flag. (default: ${PRINT})
        \\  -n <nthread> Number of threads to use. (default: ${N})
        \\  -r <right>   Right multiplicant.
        \\  -h           Display this help message and exit.
        \\EOF
        \\
        \\}; N=1; OUTPUT=null; PRINT=false;
        \\
        \\while getopts "l:n:o:pr:h" OPT; do case "$OPT" in
        \\  l ) LEFT="$OPTARG" ;;
        \\  n ) N="$OPTARG" ;;
        \\  o ) OUTPUT="$OPTARG" ;;
        \\  p ) PRINT=true ;;
        \\  r ) RIGHT="$OPTARG" ;;
        \\  h ) usage && exit 0 ;;
        \\  \?) usage && exit 1 ;;
        \\esac done
        \\
        \\[[ "$OUTPUT" != "null" ]] && OUTPUT="[\"$OUTPUT\"]"
        \\
        \\echo '{
        \\  "matrix" : {
        \\    "multiply" : {},
        \\    "inputs" : ["'"$LEFT"'", "'"$RIGHT"'"],
        \\    "outputs" : '"$OUTPUT"',
        \\    "print" : '$PRINT'
        \\  }
        \\}' > $INPUT_JSON
    ;

    const content_matsort =
        \\  -c <column> Reference column. (default: ${COLUMN})
        \\  -i <input>  Matrix to sort.
        \\  -o <output> Output matrix. (default: ${OUTPUT})
        \\  -p <print>  Boolean print flag. (default: ${PRINT})
        \\  -h          Display this help message and exit.
        \\EOF
        \\
        \\}; COLUMN=0; OUTPUT=null; PRINT=false
        \\
        \\while getopts "c:i:o:ph" OPT; do case "$OPT" in
        \\  c ) COLUMN="$OPTARG" ;;
        \\  i ) MATRIX="$OPTARG" ;;
        \\  o ) OUTPUT="$OPTARG" ;;
        \\  p ) PRINT=true ;;
        \\  h ) usage && exit 0 ;;
        \\  \?) usage && exit 1 ;;
        \\esac done
        \\
        \\[[ "$OUTPUT" != "null" ]] && OUTPUT="\"$OUTPUT\""
        \\
        \\echo '{
        \\  "sort" : {
        \\    "input" : "'"$MATRIX"'",
        \\    "algorithm" : "bubble",
        \\    "output" : '"$OUTPUT"',
        \\    "column" : '$COLUMN',
        \\    "print" : '$PRINT'
        \\  }
        \\}' > $INPUT_JSON
    ;

    const content_mersenne =
        \\  -c <count>  Number of Mersenne primes to generate. (default: ${COUNT})
        \\  -o <output> Output file. (default: ${OUTPUT})
        \\  -s <start>  Starting number. (default: ${START})
        \\  -h          Display this help message and exit.
        \\EOF
        \\
        \\}; COUNT=10; OUTPUT=null; START=1
        \\
        \\while getopts "c:o:s:h" OPT; do case "$OPT" in
        \\  c ) COUNT="$OPTARG" ;;
        \\  o ) OUTPUT="$OPTARG" ;;
        \\  s ) START="$OPTARG" ;;
        \\  h ) usage && exit 0 ;;
        \\  \?) usage && exit 1 ;;
        \\esac done
        \\
        \\[[ "$OUTPUT" != "null" ]] && OUTPUT="\"$OUTPUT\""
        \\
        \\echo '{
        \\  "prime" : {
        \\    "mode" : "mersenne",
        \\    "generate" : {
        \\      "count" : '"$COUNT"',
        \\      "output" : '"$OUTPUT"'
        \\    },
        \\    "number" : '"$START"'
        \\  }
        \\}' > $INPUT_JSON
    ;

    const content_prime =
        \\  -c <count>        Number of primes to generate. (default: ${COUNT})
        \\  -l <log_interval> Log interval for output. (default: ${LOG_INTERVAL})
        \\  -o <output>       Output file. (default: ${OUTPUT})
        \\  -s <start>        Starting number. (default: ${START})
        \\  -h                Display this help message and exit.
        \\EOF
        \\
        \\}; COUNT=10; LOG_INTERVAL=1; OUTPUT=null; START=1
        \\
        \\while getopts "c:l:o:s:h" OPT; do case "$OPT" in
        \\  c ) COUNT="$OPTARG" ;;
        \\  l ) LOG_INTERVAL="$OPTARG" ;;
        \\  o ) OUTPUT="$OPTARG" ;;
        \\  s ) START="$OPTARG" ;;
        \\  h ) usage && exit 0 ;;
        \\  \?) usage && exit 1 ;;
        \\esac done
        \\
        \\[[ "$OUTPUT" != "null" ]] && OUTPUT="\"$OUTPUT\""
        \\
        \\echo '{
        \\  "prime" : {
        \\    "mode" : "basic",
        \\    "generate" : {
        \\      "count" : '"$COUNT"',
        \\      "log_interval" : '"$LOG_INTERVAL"',
        \\      "output" : '"$OUTPUT"'
        \\    },
        \\    "number" : '"$START"'
        \\  }
        \\}' > $INPUT_JSON
    ;

    const content_randmat =
        \\  -c <columns>      Number of cols. (default: ${COLS})
        \\  -d <distribution> Distribution. (default: ${DISTRIBUTION})
        \\  -o <output>       Output matrix. (default: ${OUTPUT})
        \\  -p <print>        Boolean print flag. (default: ${PRINT})
        \\  -r <rows>         Number of rows. (default: ${ROWS})
        \\  -s <seed>         Random seed. (default: ${SEED})
        \\  -h                Display this help message and exit.
        \\EOF
        \\
        \\}; COLS=2; DISTRIBUTION="normal"; OUTPUT=null; PRINT=false; ROWS=2; SEED=1
        \\
        \\while getopts "c:d:o:pr:s:h" OPT; do case "$OPT" in
        \\  c ) COLS="$OPTARG" ;;
        \\  d ) DISTRIBUTION="$OPTARG" ;;
        \\  o ) OUTPUT="$OPTARG" ;;
        \\  p ) PRINT=true ;;
        \\  r ) ROWS="$OPTARG" ;;
        \\  s ) SEED="$OPTARG" ;;
        \\  h ) usage && exit 0 ;;
        \\  \?) usage && exit 1 ;;
        \\esac done
        \\
        \\[[ "$OUTPUT" != "null" ]] && OUTPUT="[\"$OUTPUT\"]"
        \\
        \\echo '{
        \\  "matrix" : {
        \\    "random" : {
        \\      "distribution" : "'"$DISTRIBUTION"'",
        \\      "parameters" : [0, 1],
        \\      "seed" : '$SEED',
        \\      "dims" : ['$ROWS', '$COLS']
        \\    },
        \\    "outputs" : '"$OUTPUT"',
        \\    "print" : '$PRINT'
        \\  }
        \\}' > $INPUT_JSON
    ;

    const content_schlitter =
        \\  -i <input>        Trajectory to use. (default: ${INPUT})
        \\  -t <temperature>  Temperature. (default: ${TEMP})
        \\  -h                Display this help message and exit.
        \\EOF
        \\
        \\}; INPUT="trajectory.xyz"; TEMP=298.15
        \\
        \\while getopts "i:t:h" OPT; do case "$OPT" in
        \\  i ) INPUT="$OPTARG" ;;
        \\  t ) TEMP="$OPTARG" ;;
        \\  h ) usage && exit 0 ;;
        \\  \?) usage && exit 1 ;;
        \\esac done
        \\
        \\echo '{
        \\  "trajectory_analysis" : {
        \\    "input" : "'"$INPUT"'",
        \\    "remove_translation" : true,
        \\    "remove_rotation" : true,
        \\    "entropy" : {
        \\      "algorithm" : "schlitter",
        \\      "temperature" : '$TEMP'
        \\    }
        \\  }
        \\}' > $INPUT_JSON
    ;

    _ = try file_fibonacci.write(header ++ "\n" ++ content_fibonacci ++ "\n\nexport OMP_NUM_THREADS=1; acorn $INPUT_JSON"); file_fibonacci.close();
    _ = try file_matsort  .write(header ++ "\n" ++ content_matsort   ++ "\n\nexport OMP_NUM_THREADS=1; acorn $INPUT_JSON");  file_matsort .close();
    _ = try file_mersenne .write(header ++ "\n" ++ content_mersenne  ++ "\n\nexport OMP_NUM_THREADS=1; acorn $INPUT_JSON");  file_mersenne.close();
    _ = try file_prime    .write(header ++ "\n" ++ content_prime     ++ "\n\nexport OMP_NUM_THREADS=1; acorn $INPUT_JSON");     file_prime.close();
    _ = try file_randmat  .write(header ++ "\n" ++ content_randmat   ++ "\n\nexport OMP_NUM_THREADS=1; acorn $INPUT_JSON");   file_randmat.close();
    _ = try file_schlitter.write(header ++ "\n" ++ content_schlitter ++ "\n\nexport OMP_NUM_THREADS=1; acorn $INPUT_JSON"); file_schlitter.close();

    _ = try file_hf      .write(header ++ "\n" ++ content_hf       ++ "\n\nexport OMP_NUM_THREADS=$N; acorn $INPUT_JSON"); file_hf      .close();
    _ = try file_integral.write(header ++ "\n" ++ content_integral ++ "\n\nexport OMP_NUM_THREADS=$N; acorn $INPUT_JSON"); file_integral.close();
    _ = try file_matmul  .write(header ++ "\n" ++ content_matmul   ++ "\n\nexport OMP_NUM_THREADS=$N; acorn $INPUT_JSON"); file_matmul  .close();
}
