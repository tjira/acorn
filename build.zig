const std = @import("std"); const builtin = @import("builtin");

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

        const main_executable = builder.addExecutable(.{
            .name = "acorn",
            .optimize = optimize,
            .root_source_file = builder.path("src/main.zig"),
            .single_threaded = true,
            .strip = if (optimize != .Debug) true else false,
            .target = builder.resolveTargetQuery(target)
        });

        const output = builder.addInstallArtifact(main_executable, .{
            .dest_dir = .{.override = .{.custom = try target.zigTriple(builder.allocator)}}
        });

        builder.getInstallStep().dependOn(&output.step);

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
