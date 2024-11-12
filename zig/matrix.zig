const std = @import("std");

const Vector = @import("vector.zig").Vector;

pub fn Matrix(comptime T: type) type {
    return struct {
        data: std.ArrayList(T),
        rows: usize, cols: usize,

        pub fn add(self: Matrix(T), other: Matrix(T)) !Matrix(T) {
            if (self.rows != other.rows or self.cols != other.cols) return error.IncompatibleMatrices; // throw an error if the dimensions do not match
            var result = try self.clone(); // initialize the resulting matrix with self
            for (0..self.data.items.len) |i| result.data.items[i] += other.data.items[i]; // add the other matrix to the first one
            return result; // return the resulting matrix
        }
        pub fn addScalar(self: Matrix(T), scalar: T) Matrix(T) {
            var result = Matrix(T){.data = try self.data.clone(), .rows = self.rows, .cols = self.cols}; // initialize the resulting matrix with self
            for (0..self.data.items.len) |i| result.data.items[i] += scalar; // add the scalar to the matrix
            return result; // return the resulting matrix
        }
        pub fn at(self: Matrix(T), i: usize, j: usize) !T {
            if (i >= self.rows or j >= self.cols) return error.IndexOutOfRange; // return error if the index is out of range
            return self.data.items[i * self.cols + j]; // access the correct element in the 1D list
        }
        pub fn debug(self: Matrix(T)) void {
            for (0..self.rows) |i| for (0..self.cols) |j| std.debug.print("{any}{s}", .{self.at(i, j), if (j < self.cols - 1) " " else "\n"}); // print each element of the matrix
        }
        pub fn eqApprox(self: Matrix(T), other: Matrix(T), tolerance: f64) !bool {
            if (self.rows != other.rows or self.cols != other.cols) return error.IncompatibleMatrices; // throw an error if the dimensions do not match
            for (0..self.rows) |i| for (0..self.cols) |j| if (try self.at(i, j) - try other.at(i, j) > tolerance) return false; // return false if the difference of any of the elements is above tolerance
            return true; // return true if all is below tolerance
        }
        pub fn init(m: usize, n: usize, data: []const T, allocator: std.mem.Allocator) !Matrix(T) {
            var matrix = Matrix(T){.data = std.ArrayList(T).init(allocator), .rows = m, .cols = n}; // initialize the matrix
            try matrix.data.appendSlice(data); // append the data
            return matrix; // return the resulting matrix
        }
        pub fn row(self: Matrix(T), i: usize) !Vector(T) {
            if (i >= self.rows) return error.IndexOutOfRange; // return error if the index is out of range
            var vector = try Vector(T).zero(self.rows, self.data.allocator); // initialize the output vector as zeros
            for (0..self.cols) |j| try vector.set(j, try self.at(i, j)); // set the values to the output vector
            return vector; // return the output vector
        }
        pub fn set(self: Matrix(T), i: usize, j: usize, value: T) !void {
            if (i >= self.rows or j >= self.cols) return error.IndexOutOfRange; // return error if the index is out of range
            self.data.items[i * self.cols + j] = value; // access the correct element in the 1D list and set the value
        }
        pub fn subScalar(self: Matrix(T), scalar: T) Matrix(T) {
            return self.addScalar(-scalar); // return the resulting matrix as adding the negative scalar
        }
        pub fn zero(m: usize, n: usize, allocator: std.mem.Allocator) !Matrix(T) {
            var matrix = Matrix(T){.data = std.ArrayList(T).init(allocator), .rows = m, .cols = n}; // initialize the matrix
            try matrix.data.appendNTimes(0, m * n); // append the zeros
            return matrix; // return the resulting matrix
        }
    };
}
