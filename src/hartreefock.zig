const std = @import("std");

const mat = @import("matrix.zig");

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

pub fn HartreeFockOptions(comptime T: type) type {
    return struct {
        const Integral = struct {
            overlap: []const u8 = "S_AO.mat",
            kinetic: []const u8 = "T_AO.mat",
            nuclear: []const u8 = "V_AO.mat",
            coulomb: []const u8 = "J_AO.mat"
        };

        threshold: T = 1e-8,
        maxiter: u32 = 100,

        molecule: []const u8 = "molecule.xyz", integral: Integral = .{}
    };
}

pub fn HartreeFockOutput(comptime T: type) type {
    return struct {
        C: Matrix(T), eps: Vector(T),

        pub fn init(nbf: usize, allocator: std.mem.Allocator) !HartreeFockOutput(T) {
            return HartreeFockOutput(T){
                .C = try Matrix(T).init(nbf, nbf, allocator), .eps = try Vector(T).init(nbf, allocator)
            };
        }
        pub fn deinit(self: HartreeFockOutput(T)) void {
            self.C.deinit(); self.eps.deinit();
        }
    };
}

pub fn run(comptime T: type, opt: HartreeFockOptions(T), print: bool, allocator: std.mem.Allocator) !HartreeFockOutput(T) {
    const S_AO = try mat.read(T, opt.integral.overlap, allocator);
    // const T_AO = try mat.read(T, opt.integral.kinetic, allocator);
    // const V_AO = try mat.read(T, opt.integral.nuclear, allocator);
    // const J_AO = try ten.read(T, opt.integral.coulomb, allocator);

    const nbf = S_AO.cols; const output = try HartreeFockOutput(T).init(nbf, allocator);

    // _ = T_AO;
    // _ = V_AO;
    _ = print;

    return output;
}
