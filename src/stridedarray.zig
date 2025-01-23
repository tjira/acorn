//! A strided array is a view into an array with a non-1 stride.

const std = @import("std");

/// Strided array class.
pub fn StridedArray(comptime T: type) type {
    return struct {
        data: []T, len: usize, stride: usize, zero: usize,

        /// Get the element at the specified index as a value.
        pub fn at(self: StridedArray(T), i: usize) T {
            return self.ptr(i).*;
        }

        /// Return the element at the specified index as a mutable reference.
        pub fn ptr(self: StridedArray(T), i: usize) *T {
            return &self.data[self.zero + i * self.stride];
        }
    };
}
