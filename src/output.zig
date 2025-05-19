//! Module where all the output structs are stored.

const std = @import("std");

const Vector = @import("vector.zig").Vector;
const Matrix = @import("matrix.zig").Matrix;
const System = @import("system.zig").System;
const Tensor = @import("tensor.zig").Tensor;

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
        hf: HartreeFockOutput(T), E: T, G: ?Matrix(T), H: ?Matrix(T), f: ?Vector(T),

        /// Free the memory allocated for the CI output.
        pub fn deinit(self: ConfigurationInteractionOutput(T)) void {
            self.hf.deinit();

            if (self.G != null) self.G.?.deinit();
            if (self.H != null) self.H.?.deinit();
            if (self.f != null) self.f.?.deinit();
        }
    };
}

/// The Hartree-Fock output.
pub fn HartreeFockOutput(comptime T: type) type {
    return struct {
        S_AS: Matrix(T), T_AS: Matrix(T), V_AS: Matrix(T), J_AS: Tensor(T),
        C_AS: Matrix(T), D_AS: Matrix(T), E_MS: Matrix(T), F_AS: Matrix(T),

        mulliken: ?Vector(T),

        basis: std.ArrayList(T), E: T, G: ?Matrix(T), H: ?Matrix(T), f: ?Vector(T),

        /// Free the memory allocated for the Hartree-Fock output.
        pub fn deinit(self: HartreeFockOutput(T)) void {
            self.S_AS.deinit(); self.T_AS.deinit(); self.V_AS.deinit(); self.J_AS.deinit();
            self.C_AS.deinit(); self.D_AS.deinit(); self.E_MS.deinit(); self.F_AS.deinit();

            if (self.mulliken != null) self.mulliken.?.deinit();

            self.basis.deinit();

            if (self.G != null) self.G.?.deinit();
            if (self.H != null) self.H.?.deinit();
            if (self.f != null) self.f.?.deinit();
        }
    };
}

/// The Moller-Plesset output.
pub fn MollerPlessetOutput(comptime T: type) type {
    return struct {
        hf: HartreeFockOutput(T), E: T, G: ?Matrix(T), H: ?Matrix(T), f: ?Vector(T),

        /// Free the memory allocated for the Moller-Plesset output.
        pub fn deinit(self: MollerPlessetOutput(T)) void {
            self.hf.deinit();

            if (self.G != null) self.G.?.deinit();
            if (self.H != null) self.H.?.deinit();
            if (self.f != null) self.f.?.deinit();
        }
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
