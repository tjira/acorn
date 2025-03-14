const std = @import("std"); const builtin = @import("builtin");

const examples = [_][]const u8{
    "mash_test",
    "ho_excited_states",
    "uracil_lvc",
    "surface_hopping_compare",
    "model_nonadiabatic_quantum_dynamics",
};

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
    const build_examples = builder.option(bool, "BUILD_EXAMPLES", "Build examples in the research folder") orelse false;
    const debug          = builder.option(bool, "DEBUG",          "Build everything in the debug mode"   ) orelse false;

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
                .install_dir = .prefix, .install_subdir = "../docs/code", .source_dir = main_executable.getEmittedDocs()
            });

            builder.step("docs", "Compile code documentation" ).dependOn(&docs_install                           .step);
            builder.step("run",  "Run the compiled executable").dependOn(&builder.addRunArtifact(main_executable).step);
        }

        if (build_examples) for (examples) |example| {

            const example_source_file   = try std.mem.concat(builder.allocator, u8, &[_][]const u8{"research/", example, ".zig"                        });
            const example_output_folder = try std.mem.concat(builder.allocator, u8, &[_][]const u8{try target.zigTriple(builder.allocator), "/research"});

            const example_executable = builder.addExecutable(.{
                .name = example,
                .optimize = if (debug) .Debug else .ReleaseFast,
                .root_source_file = builder.path(example_source_file),
                .single_threaded = true,
                .strip = !debug,
                .target = builder.resolveTargetQuery(target)
            });

            example_executable.root_module.addImport("acorn", builder.addModule("acorn", .{.root_source_file = builder.path("src/main.zig")}));

            const example_install = builder.addInstallArtifact(example_executable,
                .{.dest_dir = .{.override = .{.custom = example_output_folder}}
            });

            builder.getInstallStep().dependOn(&example_install.step);
        };
    }

    const test_executable = builder.addTest(.{
        .optimize = if (debug) .Debug else .ReleaseFast,
        .root_source_file = builder.path("test/main.zig"),
        .single_threaded = true,
        .strip = !debug
    });

    test_executable.root_module.addImport("acorn", builder.addModule("acorn", .{.root_source_file = builder.path("src/main.zig")}));

    builder.step("test", "Run unit tests").dependOn(&builder.addRunArtifact(test_executable).step);
}
