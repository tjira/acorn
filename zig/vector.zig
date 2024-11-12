const std = @import("std");

pub fn Vector(comptime T: type) type {
    return struct {
        data: std.ArrayList(T),
        rows: usize,

        pub fn add(self: Vector(T), other: Vector(T)) !Vector(T) {
            if (self.rows != other.rows) return error.IncompatibleVectors; // throw an error if the dimensions do not match
            var result = try self.clone(); // initialize the resulting vector with self
            for (0..self.data.items.len) |i| result.data.items[i] += other.data.items[i]; // add the other vector to the first one
            return result; // return the resulting vector
        }
        pub fn addScalar(self: Vector(T), scalar: T) !Vector(T) {
            var result = Vector(T){.data = try self.data.clone(), .rows = self.rows}; // initialize the resulting vector with self
            for (0..self.data.items.len) |i| result.data.items[i] += scalar; // add the scalar to the vector
            return result; // return the resulting vector
        }
        pub fn at(self: Vector(T), i: usize) !T {
            if (i >= self.rows) return error.IndexOutOfRange; // return error if the index is out of range
            return self.data.items[i]; // access the correct element in the 1D list
        }
        pub fn clone(self: Vector(T)) !Vector(T) {
            return Vector(T){.data = try self.data.clone(), .rows = self.rows}; // return the cloned vector
        }
        pub fn debug(self: Vector(T)) void {
            for (0..self.rows) |i| std.debug.print("{any}\n", .{self.at(i)}); // print each element of the vector
        }
        pub fn divScalar(self: Vector(T), scalar: T) !Vector(T) {
            return self.mulScalar(1 / scalar); // return the resulting vector as multiplying with the inverse
        }
        pub fn init(m: usize, data: []const T, allocator: std.mem.Allocator) !Vector(T) {
            var vector = Vector(T){.data = std.ArrayList(T).init(allocator), .rows = m}; // initialize the vector
            try vector.data.appendSlice(data); // append the data
            return vector; // return the resulting vector
        }
        pub fn mulScalar(self: Vector(T), scalar: T) !Vector(T) {
            var result = Vector(T){.data = try self.data.clone(), .rows = self.rows}; // initialize the resulting vector with self
            for (0..self.data.items.len) |i| result.data.items[i] *= scalar; // multiply the vector with the scalar
            return result; // return the resulting vector
        }
        pub fn negate(self: Vector(T)) !Vector(T) {
            var result = Vector(T){.data = try self.data.clone(), .rows = self.rows}; // initialize the resulting vector with self
            for (0..self.data.items.len) |i| result.data.items[i] *= -1; // negate the values
            return result; // return the resulting vector
        }
        pub fn set(self: Vector(T), i: usize, value: T) !void {
            if (i >= self.rows) return error.IndexOutOfRange; // return error if the index is out of range
            self.data.items[i] = value; // access the correct element in the 1D list and set the value
        }
        pub fn subScalar(self: Vector(T), scalar: T) !Vector(T) {
            return self.addScalar(-scalar); // return the resulting vector as adding the negative scalar
        }
        pub fn zero(m: usize, allocator: std.mem.Allocator) !Vector(T) {
            var vector = Vector(T){.data = std.ArrayList(T).init(allocator), .rows = m}; // initialize the vector
            try vector.data.appendNTimes(0, m); // append the zeros
            return vector; // return the resulting vector
        }
    };
}
