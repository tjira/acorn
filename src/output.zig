//! Module where all the output structs are stored.

const std = @import("std");

const Vector = @import("vector.zig").Vector;
const Matrix = @import("matrix.zig").Matrix;

/// The classical dynamics output.
pub fn ClassicalDynamicsOutput(comptime T: type) type {
    return struct {
        pop: Vector(T),

        r: Vector(T),
        p: Vector(T),

        Ekin: T,
        Epot: T,

        /// Initialize the classical dynamics output.
        pub fn init(ndim: usize, nstate: usize, allocator: std.mem.Allocator) !ClassicalDynamicsOutput(T) {
            return ClassicalDynamicsOutput(T){
                .pop = try Vector(T).init(nstate, allocator),

                .r = try Vector(T).init(ndim, allocator),
                .p = try Vector(T).init(ndim, allocator),

                .Ekin = 0,
                .Epot = 0
            };
        }

        /// Free the memory allocated for the classical dynamics output.
        pub fn deinit(self: ClassicalDynamicsOutput(T)) void {
            self.pop.deinit(); self.r.deinit(); self.p.deinit();
        }
    };
}

/// The CI output.
pub fn ConfigurationInteractionOutput(comptime T: type) type {
    return struct {
        E: T
    };
}

/// The Hartree-Fock output.
pub fn HartreeFockOutput(comptime T: type) type {
    return struct {
        C_MO: Matrix(T), D_MO: Matrix(T), E_MO: Matrix(T), F_AO: Matrix(T), E: T, VNN: T, nbf: usize, nocc: usize,

        /// Free the memory allocated for the Hartree-Fock output.
        pub fn deinit(self: HartreeFockOutput(T)) void {
            self.C_MO.deinit(); self.D_MO.deinit(); self.E_MO.deinit(); self.F_AO.deinit();
        }
    };
}

/// The Moller-Plesset output.
pub fn MollerPlessetOutput(comptime T: type) type {
    return struct {
        E: T
    };
}

/// The quantum dynamics output struct.
pub fn QuantumDynamicsOutput(comptime T: type) type {
    return struct {
        P: []Matrix(T),
        r: []Vector(T),
        p: []Vector(T),

        Ekin: []T,
        Epot: []T,

        allocator: std.mem.Allocator,

        /// Initialize the quantum dynamics output struct.
        pub fn init(ndim: usize, nstate: usize, propagations: usize, allocator: std.mem.Allocator) !QuantumDynamicsOutput(T) {
            var output = QuantumDynamicsOutput(T){
                .P = try allocator.alloc(Matrix(T), propagations),
                .r = try allocator.alloc(Vector(T), propagations),
                .p = try allocator.alloc(Vector(T), propagations),

                .Ekin = try allocator.alloc(T, propagations),
                .Epot = try allocator.alloc(T, propagations),

                .allocator = allocator,
            };

            for (0..propagations) |i| {
                output.P[i] = try Matrix(T).init(nstate, nstate, allocator);
                output.r[i] = try Vector(T).init(ndim,           allocator);
                output.p[i] = try Vector(T).init(ndim,           allocator);
            }

            return output;
        }

        /// Free the memory allocated for the quantum dynamics output struct.
        pub fn deinit(self: QuantumDynamicsOutput(T)) void {
            self.allocator.free(self.Ekin); self.allocator.free(self.Epot);

            for (0..self.P.len) |i| {
                self.P[i].deinit(); self.r[i].deinit(); self.p[i].deinit();
            }

            self.allocator.free(self.P); self.allocator.free(self.r); self.allocator.free(self.p);
        }
    };
}
