//! A vector is a one-dimensional array of elements. This module provides a generic implementation of a vector.

const std = @import("std"); const Complex = std.math.Complex;

const StridedArray = @import("stridedarray.zig").StridedArray;
const Matrix       = @import("matrix.zig").Matrix            ;

const asfloat = @import("helper.zig").asfloat;
const istruct = @import("helper.zig").istruct;

/// Vector class.
pub fn Vector(comptime T: type) type {
    return struct {
        data: []T, rows: usize, allocator: std.mem.Allocator,

        /// Initialize a new column vector with the specified number of rows. Returns an error if the allocation fails.
        pub fn init(rows: usize, allocator: std.mem.Allocator) !Vector(T) {
            return Vector(T){
                .data = try allocator.alloc(T, rows),
                .rows = rows,
                .allocator = allocator
            };
        }

        /// Free the memory allocated for the vector.
        pub fn deinit(self: Vector(T)) void {
            self.allocator.free(self.data);
        }

        /// Get the element at the specified index as a value.
        pub fn at(self: Vector(T), i: usize) T {
            return self.ptr(i).*;
        }

        /// Clone the vector. Returns an error if the allocation fails.
        pub fn clone(self: Vector(T)) !Vector(T) {
            const other = try Vector(T).init(self.rows, self.allocator);

            @memcpy(other.data, self.data);

            return other;
        }

        /// Convert the vector to a complex vector. Returns an error if the allocation fails.
        pub fn complex(self: Vector(T)) !Vector(Complex(T)) {
            var other = try Vector(Complex(T)).init(self.rows, self.allocator);

            for (0..self.data.len) |i| {
                other.data[i] = Complex(T).init(self.data[i], 0);
            }

            return other;
        }

        /// Fill the vector with the specified value.
        pub fn fill(self: Vector(T), value: T) void {
            for (0..self.data.len) |i| {
                self.data[i] = value;
            }
        }

        /// Return the vector as a matrix with one column. No memory is allocated.
        pub fn matrix(self: Vector(T)) Matrix(T) {
            return Matrix(T){
                .data = self.data,
                .rows = self.rows,
                .cols = 1,
                .allocator = self.allocator
            };
        }

        /// Return the element at the specified index as a mutable reference.
        pub fn ptr(self: Vector(T), i: usize) *T {
            return &self.data[i];
        }

        /// Fill the vector with normally distributed random numbers.
        pub fn randn(self: Vector(T), mean: T, stdev: T, seed: usize) void {
            var rng = std.rand.DefaultPrng.init(seed); const random = rng.random();

            for (0..self.rows) |i| {
                self.data[i] = mean + stdev * random.floatNorm(T);
            }
        }
    };
}

/// Add two vectors. The output vector is stored in the vector w.
pub fn add(comptime T: type, w: *Vector(T), u: Vector(T), v: Vector(T)) void {
    if (comptime !istruct(T)) for (0..u.rows) |i| {
        w.ptr(i).* = u.at(i) + v.at(i);
    };

    if (comptime istruct(T)) for (0..u.rows) |i| {
        w.ptr(i).* = u.at(i).add(v.at(i));
    };
}

/// Multiply two vectors element-wise. The output vector is stored in the vector w.
pub fn mul(comptime T: type, w: *Vector(T), u: Vector(T), v: Vector(T)) void {
    if (comptime !istruct(T)) for (0..u.rows) |i| {
        w.ptr(i).* = u.at(i) * v.at(i);
    };

    if (comptime istruct(T)) for (0..u.rows) |i| {
        w.ptr(i).* = u.at(i).mul(v.at(i));
    };
}

/// Multiply a vector by a scalar. The output vector is stored in the vector w.
pub fn muls(comptime T: type, w: *Vector(T), u: Vector(T), s: T) void {
    if (comptime !istruct(T)) for (0..u.rows) |i| {
        w.ptr(i).* = u.at(i) * s;
    };

    if (comptime istruct(T)) for (0..u.rows) |i| {
        w.ptr(i).* = u.at(i).add(s);
    };
}
