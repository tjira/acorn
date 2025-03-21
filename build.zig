const std = @import("std"); const builtin = @import("builtin");

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
    const debug = builder.option(bool, "DEBUG", "Build everything in the debug mode") orelse false;

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
                .install_dir = .prefix, .install_subdir = "docs", .source_dir = builder.addExecutable(.{
                    .name = "main", .root_source_file = builder.path("src/main.zig"), .target = builder.host
                }).getEmittedDocs(),
            });

            var script_step = builder.step("script", "Generate executable scripts"); script_step.makeFn = script; script_step.dependOn(builder.getInstallStep());

            builder.step("docs", "Compile code documentation" ).dependOn(&docs_install                           .step);
            builder.step("run",  "Run the compiled executable").dependOn(&builder.addRunArtifact(main_executable).step);
        }
    }

    const test_executable = builder.addTest(.{
        .name = "test",
        .optimize = if (debug) .Debug else .ReleaseFast,
        .root_source_file = builder.path("test/main.zig"),
        .single_threaded = true,
        .strip = !debug,
        .target = builder.host
    });

    test_executable.root_module.addImport("acorn", builder.addModule("acorn", .{.root_source_file = builder.path("src/main.zig")}));

    builder.step("test", "Run unit tests").dependOn(&builder.addRunArtifact(test_executable).step);
}

pub fn script(self: *std.Build.Step, progress: std.Progress.Node) !void {
    _ = self; _ = progress;

    for (targets) |target| if (target.os_tag == .linux) {
        try linuxScripts(std.heap.page_allocator, try target.zigTriple(std.heap.page_allocator));
    };
}

pub fn linuxScripts(allocator: std.mem.Allocator, target: []const u8) !void {
    const path = try std.mem.concat(allocator, u8, &[_][]const u8{"zig-out/", target, "/"});

    const file_matsort  = try std.fs.cwd().createFile(try std.mem.concat(allocator, u8, &[_][]const u8{path, "matsort" }), .{.mode=0o755});
    const file_mersenne = try std.fs.cwd().createFile(try std.mem.concat(allocator, u8, &[_][]const u8{path, "mersenne"}), .{.mode=0o755});

    const content_matsort =
        \\#!/usr/bin/bash
        \\
        \\if [ "$#" -ne 1 ]; then
        \\  echo "USAGE: matsort MATRIX"; exit 1
        \\fi
        \\
        \\clean() {{ rm -f input.json; }}; trap clean SIGINT
        \\
        \\echo '{{
        \\  "sort" : {{
        \\    "input" : "'"$1"'",
        \\    "algorithm" : "bubble",
        \\    "output" : "'"$1.sorted"'"
        \\  }}
        \\}}' > input.json
        \\
        \\acorn
        \\
        \\mv "$1.sorted" "$1" && clean
    ;

    const content_mersenne =
        \\#!/usr/bin/bash
        \\
        \\if [ "$#" -ne 1 ]; then
        \\  echo "USAGE: mersenne COUNT"; exit 1
        \\fi
        \\
        \\clean() {{ rm -f input.json; }}; trap clean SIGINT
        \\
        \\echo '{{
        \\  "prime" : {{
        \\    "mode" : "mersenne",
        \\    "generate" : {{
        \\      "count" : '"$1"'
        \\    }},
        \\    "number" : 1
        \\  }}
        \\}}' > input.json
        \\
        \\acorn
        \\
        \\clean
    ;

    try file_matsort .writer().print(content_matsort,  .{}); file_matsort .close();
    try file_mersenne.writer().print(content_mersenne, .{}); file_mersenne.close();
}
