const std = @import("std");

pub fn StridedArray(comptime T: type) type {
    return struct {
        data: []T, len: usize, stride: usize, zero: usize,

        pub fn at(self: StridedArray(T), i: usize) T {
            return self.ptr(i).*;
        }
        pub fn ptr(self: StridedArray(T), i: usize) *T {
            return &self.data[self.zero + i * self.stride];
        }
    };
}
