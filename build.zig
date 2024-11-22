const std = @import("std");

pub fn build(builder: *std.Build) void {
    const optimize = builder.standardOptimizeOption(.{});
    const target   = builder.standardTargetOptions (.{});

    const executable = builder.addExecutable(.{
        .name = "acorn",
        .optimize = optimize,
        .root_source_file = builder.path("zig/main.zig"),
        .single_threaded = true,
        .strip = if (optimize != .Debug) true else false,
        .target = target
    });

    executable.linkSystemLibrary2("gsl", .{}); executable.linkLibC();

    builder.installArtifact(executable);
}
