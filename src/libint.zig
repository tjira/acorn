//! This file contains wrappers for the libint2 library.

const libint = @cImport(@cInclude("libint.h"));

const System = @import("system.zig").System;

/// Calculate the Coulomb integrals.
pub fn coulomb(ints: []f64, system: System(f64), basis: []const f64) void {
    libint.coulomb(&ints[0], @intCast(system.atoms.rows), &system.atoms.data[0], &system.coords.data[0], @intCast(basis.len), &basis[0]);
}

/// Calculate the kinetic integrals.
pub fn kinetic(ints: []f64, system: System(f64), basis: []const f64) void {
    libint.kinetic(&ints[0], @intCast(system.atoms.rows), &system.atoms.data[0], &system.coords.data[0], @intCast(basis.len), &basis[0]);
}

/// Calculate the nuclear integrals.
pub fn nuclear(ints: []f64, system: System(f64), basis: []const f64) void {
    libint.nuclear(&ints[0], @intCast(system.atoms.rows), &system.atoms.data[0], &system.coords.data[0], @intCast(basis.len), &basis[0]);
}

/// Calculate the overlap integrals.
pub fn overlap(ints: []f64, system: System(f64), basis: []const f64) void {
    libint.overlap(&ints[0], @intCast(system.atoms.rows), &system.atoms.data[0], &system.coords.data[0], @intCast(basis.len), &basis[0]);
}
