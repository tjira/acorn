const std = @import("std");

const hfm = @import("hartreefock.zig");

pub fn MollerPlessetOptions(comptime T: type) type {
    return struct {
        order: u32 = 2,

        hartree_fock: hfm.HartreeFockOptions(T) = .{},
    };
}

pub fn MollerPlessetOutput(comptime T: type) type {
    return struct {
        E: T,

        pub fn deinit(self: MollerPlessetOutput(T)) void {
            _ = self;
        }
    };
}

pub fn run(comptime T: type, opt: MollerPlessetOptions(T), print: bool, allocator: std.mem.Allocator) !MollerPlessetOutput(T) {
    const hf = try hfm.run(T, opt.hartree_fock, print, allocator);

    return .{
        .E = hf.E,
    };
}
