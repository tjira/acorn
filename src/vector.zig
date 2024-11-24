const std = @import("std");

pub fn Vector(comptime T: type) type {
    return struct {
        data: []T, rows: usize, allocator: std.mem.Allocator,

        pub fn init(rows: usize, allocator: std.mem.Allocator) !Vector(T) {
            return Vector(T){.data = try allocator.alloc(T, rows), .rows = rows, .allocator = allocator};
        }
        pub fn deinit(self: Vector(T)) void {
            self.allocator.free(self.data);
        }

        pub fn at(self: Vector(T), i: usize) T {
            return self.data[i];
        }
        pub fn ptr(self: Vector(T), i: usize) *T {
            return &self.data[i];
        }

        pub fn fill(self: Vector(T), value: T) void {
            for (0..self.data.len) |i| self.data[i] = value;
        }
    };
}
