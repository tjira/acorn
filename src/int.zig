//! This file contains wrappers for the libint2 library.

const int = @cImport(@cInclude("int.h"));

const std = @import("std");

/// Initialize the libint2 library.
pub fn initialize() void {
    int.initialize();
}

/// Finalize the libint2 library.
pub fn finalize() void {
    int.finalize();
}
