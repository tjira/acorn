const std = @import("std");

const hfm = @import("hartreefock.zig");

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

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

    var dets = std.ArrayList(Vector(usize)).init(allocator); defer dets.deinit();

    for (dets.items) |*e| e.deinit();

    return .{
        .E = hf.E,
    };
}
