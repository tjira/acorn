const std = @import("std");

pub fn build(b: *std.Build) void {
    const optimize = b.standardOptimizeOption(.{});
    const target   = b.standardTargetOptions (.{});

    const main_exe = b.addExecutable(.{
        .name = "acorn",
        .optimize = optimize,
        .root_source_file = b.path("zig/main.zig"),
        .single_threaded = true,
        .strip = if (optimize != .Debug) true else false,
        .target = target
    });

    main_exe.linkSystemLibrary2("fftw3", .{}); 
    main_exe.linkSystemLibrary2("gsl",   .{}); 
    main_exe.linkLibC();
    // main_exe.linkLibCpp();

    // const test_exe = b.addTest(.{
    //     .root_source_file = b.path("test/main.zig"),
    //     .single_threaded = true,
    //     .strip = if (optimize != .Debug) true else false,
    //     .target = target
    // });

    b.installArtifact(main_exe);

    // b.step("test", "Run the test executable").dependOn(&b.addRunArtifact(test_exe).step);
    b.step("run",  "Run the main executable").dependOn(&b.addRunArtifact(main_exe).step);

}
