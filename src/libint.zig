//! This file contains wrappers for the libint2 library.

const libint = @cImport(@cInclude("libint.h"));

const std = @import("std");

/// Calculate the Coulomb integrals.
pub fn coulomb(ints: []f64, system: []const u8, basis: []const u8) void {
    libint.coulomb(&ints[0], &system[0], &basis[0]);
}

/// Calculate the kinetic integrals.
pub fn kinetic(ints: []f64, system: []const u8, basis: []const u8) void {
    libint.kinetic(&ints[0], &system[0], &basis[0]);
}

/// Calculate the nuclear integrals.
pub fn nuclear(ints: []f64, system: []const u8, basis: []const u8) void {
    libint.nuclear(&ints[0], &system[0], &basis[0]);
}

/// Calculate the overlap integrals.
pub fn overlap(ints: []f64, system: []const u8, basis: []const u8) void {
    libint.overlap(&ints[0], &system[0], &basis[0]);
}
