//! This module provides a tensor class and functions to manipulate tensors.

const std = @import("std");

const mth = @import("math.zig");

const Matrix = @import("matrix.zig").Matrix;

/// Tensor class.
pub fn Tensor(comptime T: type) type {
    return struct {
        data: []T, shape: []usize, stride: []usize, allocator: std.mem.Allocator,

        /// Initialize a new tensor with the specified shape. Returns an error if the allocation fails.
        pub fn init(shape: []const usize, allocator: std.mem.Allocator) !Tensor(T) {
            const ten = Tensor(T){
                .data = try allocator.alloc(T, mth.prod(usize, shape)),
                .shape = try allocator.alloc(usize, shape.len),
                .stride = try allocator.alloc(usize, shape.len),
                .allocator = allocator
            };

            @memcpy(ten.shape, shape);

            for (0..shape.len) |i| {
                ten.stride[i] = mth.prod(usize, shape[i + 1..shape.len]);
            }

            return ten;
        }

        /// Free the memory allocated for the tensor.
        pub fn deinit(self: Tensor(T)) void {
            self.allocator.free(self.data  );
            self.allocator.free(self.shape );
            self.allocator.free(self.stride);
        }

        /// Returns the element at the specified indices as a value.
        pub fn at(self: Tensor(T), indices: []const usize) T {
            return self.ptr(indices).*;
        }

        /// Fill the tensor with the specified value.
        pub fn fill(self: *Tensor(T), value: T) void {
            @memset(self.data, value);
        }

        /// Returns the tensor as a matrix. No memory is allocated. Currently only works for 4D tensors.
        pub fn matrix(self: Tensor(T)) Matrix(T) {
            return Matrix(T){
                .data = self.data,
                .rows = self.shape[0] * self.shape[1],
                .cols = self.shape[2] * self.shape[3],
                .allocator = self.allocator
            };
        }

        /// Print the matrix to the given device.
        pub fn print(self: Tensor(T), device: anytype) !void {
            var buffered = std.io.bufferedWriter(device); var writer = buffered.writer();

            for (0..self.shape.len) |i| {
                try writer.print("{d}{s}", .{self.shape[i], if(i < self.shape.len - 1) " " else "\n"});
            }

            for (self.data, 1..) |e, i| {
                try writer.print("{d:20.14}{s}", .{e, if(i % (self.shape[2] * self.shape[3]) == 0) "\n" else " "});
            }

            try buffered.flush();
        }

        /// Returns the element at the specified indices as a mutable reference.
        pub fn ptr(self: Tensor(T), indices: []const usize) *T {
            var index: usize = 0;

            for (0..indices.len) |i| {
                index += indices[i] * self.stride[i];
            }

            return &self.data[index];
        }

        /// Write the tensor to a file.
        pub fn write(self: Tensor(T), path: []const u8) !void {
            const file = try std.fs.cwd().createFile(path, .{}); defer file.close();

            try self.print(file.writer());
        }
    };
}

/// Add two tensors element-wise. The result is stored in the tensor C.
pub fn add(comptime T: type, C: *Tensor(T), A: Tensor(T), B: Tensor(T)) void {
    if (@typeInfo(T) != .Struct) for (0..C.data.len) |i| {
        C.data[i] = A.data[i] + B.data[i];
    };

    if (@typeInfo(T) == .Struct) for (0..C.data.len) |i| {
        C.data[i] = A.data[i].add(B.data[i]);
    };
}

/// Multiply two tensors element-wise. The result is stored in the tensor C.
pub fn mul(comptime T: type, C: *Tensor(T), A: Tensor(T), B: Tensor(T)) void {
    if (@typeInfo(T) != .Struct) for (0..C.data.len) |i| {
        C.data[i] = A.data[i] * B.data[i];
    };

    if (@typeInfo(T) == .Struct) for (0..C.data.len) |i| {
        C.data[i] = A.data[i].mul(B.data[i]);
    };
}

/// Multiply a tensor by a scalar. The result is stored in the tensor C.
pub fn muls(comptime T: type, C: *Tensor(T), A: Tensor(T), s: T) void {
    if (@typeInfo(T) != .Struct) for (0..C.data.len) |i| {
        C.data[i] = A.data[i] * s;
    };

    if (@typeInfo(T) == .Struct) for (0..C.data.len) |i| {
        C.data[i] = A.data[i].mul(s);
    };
}

/// Read a tensor from a file. The tensor is stored in row-major order.
pub fn read(comptime T: type, path: []const u8, dim: usize, allocator: std.mem.Allocator) !Tensor(T) {
    const file = try std.fs.cwd().openFile(path, .{}); defer file.close();

    var shape = try allocator.alloc(usize, dim); defer allocator.free(shape);

    var freader = std.io.bufferedReader(file.reader()); var reader = freader.reader();
    var hbuffer: [16]u8 = undefined; var hstream = std.io.fixedBufferStream(&hbuffer);

    for (0..dim) |i| {

        hstream.reset();

        try reader.streamUntilDelimiter(hstream.writer(), if (i < dim - 1) ' ' else '\n',  4);

        shape[i] = try std.fmt.parseInt(usize, hbuffer[0..@intCast(try hstream.getPos())], 10);
    }

    const ten = try Tensor(T).init(shape, allocator);

    for (0..mth.prod(usize, shape)) |i| {

        const bytes = try reader.readBytesNoEof(20); try reader.skipBytes(1, .{});

        var j: u32 = 0; while (bytes[j] == ' ') : (j += 1) {}

        ten.data[i] = try std.fmt.parseFloat(T, bytes[j..bytes.len]);
    }

    return ten;
}

/// Subtract two tensors element-wise. The result is stored in the tensor C.
pub fn sub(comptime T: type, C: *Tensor(T), A: Tensor(T), B: Tensor(T)) void {
    if (@typeInfo(T) != .Struct) for (0..C.data.len) |i| {
        C.data[i] = A.data[i] - B.data[i];
    };

    if (@typeInfo(T) == .Struct) for (0..C.data.len) |i| {
        C.data[i] = A.data[i].sub(B.data[i]);
    };
}

/// Transpose a tensor. The result is stored in the tensor B. The axes argument specifies the permutation of the axes. Currently only works for 4D tensors.
pub fn transpose(comptime T: type, B: *Tensor(T), A: Tensor(T), axes: []const usize) void {
    if (axes.len == 4) for (0..A.shape[0]) |i| for (0..A.shape[1]) |j| for (0..A.shape[2]) |k| for (0..A.shape[3]) |l| {

        const source = [_]usize{i, j, k, l}; var destination: [4]usize = undefined;

        for (0..axes.len) |m| destination[m] = source[axes[m]];

        B.ptr(&destination).* = A.at(&source);
    }; 
}
