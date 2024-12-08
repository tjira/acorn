const std = @import("std");

const targets: []const std.Target.Query = &.{
    .{.os_tag = .linux  , .cpu_arch = .aarch64},
    .{.os_tag = .linux  , .cpu_arch = .x86_64 },
    .{.os_tag = .windows, .cpu_arch = .aarch64},
    .{.os_tag = .windows, .cpu_arch = .x86_64 },
    .{.os_tag = .macos  , .cpu_arch = .aarch64},
    .{.os_tag = .macos  , .cpu_arch = .x86_64 },
};

pub fn build(builder: *std.Build) !void {
    const optimize = builder.standardOptimizeOption(.{});

    for (targets) |target| {

        const executable = builder.addExecutable(.{
            .name = "acorn",
            .optimize = optimize,
            .root_source_file = builder.path("src/main.zig"),
            .single_threaded = true,
            .strip = if (optimize != .Debug) true else false,
            .target = builder.resolveTargetQuery(target)
        });

        const output = builder.addInstallArtifact(executable, .{
            .dest_dir = .{.override = .{.custom = try target.zigTriple(builder.allocator)}}
        });

        builder.getInstallStep().dependOn(&output.step);
    }

    const test_filters = builder.option([]const []const u8, "test", "Select only specific tests") orelse &[0][]const u8{};

    const test_executable = builder.addTest(.{
        .filters = test_filters,
        .optimize = optimize,
        .root_source_file = builder.path("test/main.zig"),
        .single_threaded = true,
        .strip = if (optimize != .Debug) true else false
    });

    test_executable.root_module.addImport("acorn", builder.addModule("acorn", .{.root_source_file = builder.path("src/main.zig")}));

    builder.step("test", "Run unit tests").dependOn(&builder.addRunArtifact(test_executable).step);
}
