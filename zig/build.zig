const std = @import("std");

pub fn build(builder: *std.Build) void {
    builder.installArtifact(builder.addExecutable(.{.name = "acorn", .optimize = .Debug, .root_source_file = builder.path("main.zig"), .target = builder.host}));
}
