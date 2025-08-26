const std = @import("std"); const builtin = @import("builtin");

const targets: []const std.Target.Query = &.{
    .{.os_tag = .linux, .cpu_arch = .aarch64, .abi = .gnu },
    .{.os_tag = .linux, .cpu_arch = .aarch64, .abi = .musl},
    .{.os_tag = .linux, .cpu_arch = .x86_64,  .abi = .gnu },
    .{.os_tag = .linux, .cpu_arch = .x86_64,  .abi = .musl},
};

pub fn build(builder: *std.Build) !void {
    const debug = builder.option(bool, "DEBUG", "Build everything in the debug mode") orelse false;

    for (targets) |target| {

        const main_executable = builder.addExecutable(.{
            .name = "acorn",
            .root_module = builder.createModule(.{
                .optimize = if (debug) .Debug else .ReleaseFast,
                .root_source_file = builder.path("src/main.zig"),
                .strip = !debug,
                .target = builder.resolveTargetQuery(target)
            })
        });

        const main_executable_install = builder.addInstallArtifact(main_executable, .{
            .dest_dir = .{.override = .{.custom = try target.zigTriple(builder.allocator)}}
        });

        builder.getInstallStep().dependOn(&main_executable_install.step);

        if (builtin.target.cpu.arch == target.cpu_arch and builtin.target.os.tag == target.os_tag and target.abi == .musl) {

            builder.step("run", "Run the compiled executable").dependOn(&builder.addRunArtifact(main_executable).step);

            const docs = builder.addInstallDirectory(.{.source_dir = main_executable.getEmittedDocs(), .install_dir = .{.custom = "../docs"}, .install_subdir = "code"});

            builder.step("docs", "Install the code documentation").dependOn(&docs.step);
        }
    }
}
