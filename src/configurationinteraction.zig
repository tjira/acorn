const std = @import("std");

const hfm = @import("hartreefock.zig");

pub fn ConfigurationInteractionOptions(comptime T: type) type {
    return struct {
        excitation: ?[]const u32 = null,

        hartree_fock: hfm.HartreeFockOptions(T) = .{},
    };
}

pub fn ConfigurationInteractionOutput(comptime T: type) type {
    return struct {
        E: T,

        pub fn deinit(self: ConfigurationInteractionOutput(T)) void {
            _ = self;
        }
    };
}

pub fn run(comptime T: type, opt: ConfigurationInteractionOptions(T), print: bool, allocator: std.mem.Allocator) !ConfigurationInteractionOutput(T) {
    const hf = try hfm.run(T, opt.hartree_fock, print, allocator);

    return .{
        .E = hf.E,
    };
}
