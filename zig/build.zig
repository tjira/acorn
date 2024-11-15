const std = @import("std");

pub fn build(builder: *std.Build) void {
    const optimize = builder.standardOptimizeOption(.{});

    const exe = builder.addExecutable(.{
        .name = "acorn",
        .optimize = optimize,
        .root_source_file = builder.path("main.zig"),
        .target = builder.host,
        .strip = true
    });

    builder.installArtifact(exe);

    builder.step("run", "Run the application").dependOn(&builder.addRunArtifact(exe).step);
}
