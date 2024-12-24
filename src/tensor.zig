const std = @import("std");

const Matrix = @import("matrix.zig").Matrix;

const prod = @import("helper.zig").prod;

pub fn Tensor(comptime T: type) type {
    return struct {
        data: []T, shape: []usize, stride: []usize, allocator: std.mem.Allocator,

        pub fn init(shape: []const usize, allocator: std.mem.Allocator) !Tensor(T) {
            const ten = Tensor(T){
                .data = try allocator.alloc(T, prod(usize, shape)),
                .shape = try allocator.alloc(usize, shape.len),
                .stride = try allocator.alloc(usize, shape.len),
                .allocator = allocator
            };

            @memcpy(ten.shape, shape);

            for (0..shape.len) |i| {
                ten.stride[i] = prod(usize, shape[i + 1..shape.len]);
            }

            return ten;
        }
        pub fn deinit(self: Tensor(T)) void {
            self.allocator.free(self.data );
            self.allocator.free(self.shape);
        }

        pub fn at(self: Tensor(T), indices: []const usize) T {
            return self.ptr(indices).*;
        }
        pub fn ptr(self: Tensor(T), indices: []const usize) *T {
            var index: usize = 0;

            for (0..indices.len) |i| {
                index += indices[i] * self.stride[i];
            }

            return &self.data[index];
        }
        pub fn matrix(self: Tensor(T)) Matrix(T) {
            return Matrix(T){
                .data = self.data,
                .rows = self.shape[0] * self.shape[1],
                .cols = self.shape[2] * self.shape[3],
                .allocator = self.allocator
            };
        }

        pub fn fill(self: *Tensor(T), value: T) void {
            for (0..self.data.len) |i| {
                self.data[i] = value;
            }
        }
    };
}

pub fn add(comptime T: type, C: *Tensor(T), A: Tensor(T), B: Tensor(T)) void {
    if (@typeInfo(T) != .Struct) for (0..C.data.len) |i| {
        C.data[i] = A.data[i] + B.data[i];
    };

    if (@typeInfo(T) == .Struct) for (0..C.data.len) |i| {
        C.data[i] = A.data[i].add(B.data[i]);
    };
}

pub fn mul(comptime T: type, C: *Tensor(T), A: Tensor(T), B: Tensor(T)) void {
    if (@typeInfo(T) != .Struct) for (0..C.data.len) |i| {
        C.data[i] = A.data[i] * B.data[i];
    };

    if (@typeInfo(T) == .Struct) for (0..C.data.len) |i| {
        C.data[i] = A.data[i].mul(B.data[i]);
    };
}

pub fn muls(comptime T: type, C: *Tensor(T), A: Tensor(T), s: T) void {
    if (@typeInfo(T) != .Struct) for (0..C.data.len) |i| {
        C.data[i] = A.data[i] * s;
    };

    if (@typeInfo(T) == .Struct) for (0..C.data.len) |i| {
        C.data[i] = A.data[i].mul(s);
    };
}

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

    for (0..prod(usize, shape)) |i| {

        const bytes = try reader.readBytesNoEof(20); try reader.skipBytes(1, .{});

        var j: u32 = 0; while (bytes[j] == ' ') : (j += 1) {}

        ten.data[i] = try std.fmt.parseFloat(T, bytes[j..bytes.len]);
    }

    return ten;
}

pub fn sub(comptime T: type, C: *Tensor(T), A: Tensor(T), B: Tensor(T)) void {
    if (@typeInfo(T) != .Struct) for (0..C.data.len) |i| {
        C.data[i] = A.data[i] - B.data[i];
    };

    if (@typeInfo(T) == .Struct) for (0..C.data.len) |i| {
        C.data[i] = A.data[i].sub(B.data[i]);
    };
}

pub fn transpose(comptime T: type, B: *Tensor(T), A: Tensor(T), axes: []const usize) void {
    if (axes.len == 4) for (0..A.shape[0]) |i| for (0..A.shape[1]) |j| for (0..A.shape[2]) |k| for (0..A.shape[3]) |l| {

        const source = [_]usize{i, j, k, l}; var destination: [4]usize = undefined;

        for (0..axes.len) |m| destination[m] = source[axes[m]];

        B.ptr(&destination).* = A.at(&source);
    }; 
}
