//! A strided array is a view into an array with a non-1 stride.

const std = @import("std");

const asfloat = @import("helper.zig").asfloat;

/// Strided array class.
pub fn StridedArray(comptime T: type) type {
    return struct {
        data: []T, len: usize, stride: usize, zero: usize,

        /// Return the absolute sum of the strided array.
        pub fn absSum(self: StridedArray(T)) T {
            var result: T = 0;

            for (0..self.len) |i| {
                result += @abs(self.at(i));
            }

            return result;
        }

        /// Get the element at the specified index as a value.
        pub fn at(self: StridedArray(T), i: usize) T {
            return self.ptr(i).*;
        }

        /// Fill the strided array elements with liearly spaced values. Start and end are inclusive.
        pub fn linspace(self: StridedArray(T), start: T, end: T) void {
            self.ptr(0).* = start;

            for (1..self.len) |i| {
                self.ptr(i).* = start + asfloat(T, i) * (end - start) / asfloat(T, self.len - 1);
            }
        }

        /// Return the element at the specified index as a mutable reference.
        pub fn ptr(self: StridedArray(T), i: usize) *T {
            return &self.data[self.zero + i * self.stride];
        }

        /// Return the sum of the strided array.
        pub fn sum(self: StridedArray(T)) T {
            var result: T = 0;

            for (0..self.len) |i| {
                result += self.at(i);
            }

            return result;
        }
    };
}
