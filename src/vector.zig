//! A vector is a one-dimensional array of elements. This module provides a generic implementation of a vector.

const std = @import("std"); const Complex = std.math.Complex;

const StridedArray = @import("stridedarray.zig").StridedArray;
const Matrix       = @import("matrix.zig").Matrix            ;

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

        /// Get the element at the specified index as a value.
        pub fn at(self: Vector(T), i: usize) T {
            return self.ptr(i).*;
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

        /// Return the vector as a strided array with a stride of 1. No memory is allocated.
        pub fn sa(self: Vector(T)) StridedArray(T) {
            return StridedArray(T){.data = self.data, .len = self.rows, .stride = 1, .zero = 0};
        }

        /// Fill the vector with the specified value.
        pub fn fill(self: Vector(T), value: T) void {
            for (0..self.data.len) |i| {
                self.data[i] = value;
            }
        }
    };
}

/// Compare two vectors element-wise. Returns true if the vectors are equal within the specified epsilon.
pub fn eq(comptime T: type, u: Vector(T), v: Vector(T), epsilon: f64) bool {
    if (u.rows != v.rows) return false;

    if (@typeInfo(T) != .Struct) for (0..u.data.len) |i| if (@abs(u.data[i] - v.data[i]) > epsilon or std.math.isNan(@abs(u.data[i] - v.data[i]))) {
        return false;
    };

    if (@typeInfo(T) == .Struct) for (0..u.data.len) |i| if (@abs(u.data[i].re - v.data[i].re) > epsilon or @abs(u.data[i].im - v.data[i].im) > epsilon) {
        return false;
    };

    return true;
}
