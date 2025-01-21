const std = @import("std"); const builtin = @import("builtin");

const examples: [4][]const u8 = [_][]const u8{
    "example/research/fssh_vs_mash_identity_check.zig",
    "example/research/ho_excited_states.zig",
    "example/research/uracil_lvc.zig",
    "example/research/surface_hopping_compare.zig",
};

const targets: []const std.Target.Query = &.{
    .{.os_tag = .linux  , .cpu_arch = .aarch64},
    .{.os_tag = .linux  , .cpu_arch = .arm    },
    .{.os_tag = .linux  , .cpu_arch = .riscv64},
    .{.os_tag = .linux  , .cpu_arch = .x86    },
    .{.os_tag = .linux  , .cpu_arch = .x86_64 },
    .{.os_tag = .windows, .cpu_arch = .aarch64},
    .{.os_tag = .windows, .cpu_arch = .x86    },
    .{.os_tag = .windows, .cpu_arch = .x86_64 },
    .{.os_tag = .macos  , .cpu_arch = .aarch64},
    .{.os_tag = .macos  , .cpu_arch = .x86_64 },
};

pub fn build(builder: *std.Build) !void {
    const optimize = builder.standardOptimizeOption(.{});

    for (targets) |target| {

        const target_name = try target.zigTriple(builder.allocator);

        const main_executable = builder.addExecutable(.{
            .name = "acorn",
            .optimize = optimize,
            .root_source_file = builder.path("src/main.zig"),
            .single_threaded = true,
            .strip = if (optimize != .Debug) true else false,
            .target = builder.resolveTargetQuery(target)
        });

        const main_output = builder.addInstallArtifact(main_executable, .{
            .dest_dir = .{.override = .{.custom = target_name}}
        });

        builder.getInstallStep().dependOn(&main_output.step);

        for (examples) |example| {

            var bin_folder = try builder.allocator.alloc(u8, target_name.len + 8); bin_folder[target_name.len] = '/';

            @memcpy(bin_folder[0                  ..target_name.len    ], target_name      );
            @memcpy(bin_folder[target_name.len + 1..target_name.len + 8], example    [0..7]);

            const example_executable = builder.addExecutable(.{
                .name = example[17..example.len - 4],
                .optimize = optimize,
                .root_source_file = builder.path(example),
                .single_threaded = true,
                .strip = if (optimize != .Debug) true else false,
                .target = builder.resolveTargetQuery(target)
            });

            example_executable.root_module.addImport("acorn", builder.addModule("acorn", .{.root_source_file = builder.path("src/main.zig")}));

            const example_output = builder.addInstallArtifact(example_executable, .{
                .dest_dir = .{.override = .{.custom = bin_folder}}
            });

            builder.getInstallStep().dependOn(&example_output.step);
        }

        if (builtin.target.cpu.arch == target.cpu_arch and builtin.target.os.tag == target.os_tag) {

            const run_executable = builder.addRunArtifact(main_executable);

            builder.step("run",  "Run the executable").dependOn(&run_executable.step);
        }
    }

    const docs = builder.addInstallDirectory(.{
        .install_dir = .prefix,
        .install_subdir = "code",
        .source_dir = builder.addExecutable(.{
            .name = "main",
            .root_source_file = builder.path("src/main.zig"),
            .target = builder.host
        }).getEmittedDocs(),
    });

    const test_executable = builder.addTest(.{
        .optimize = optimize,
        .root_source_file = builder.path("test/main.zig"),
        .single_threaded = true,
        .strip = if (optimize != .Debug) true else false
    });

    test_executable.root_module.addImport("acorn", builder.addModule("acorn", .{.root_source_file = builder.path("src/main.zig")}));

    builder.step("docs", "Compile documentation").dependOn(&docs                                   .step);
    builder.step("test", "Run unit tests"       ).dependOn(&builder.addRunArtifact(test_executable).step);
}
