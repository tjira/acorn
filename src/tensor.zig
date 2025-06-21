//! This module provides a tensor class and functions to manipulate tensors.

const std = @import("std"); const cwp = @import("cwrapper.zig");

const mth = @import("math.zig");

const Matrix = @import("matrix.zig").Matrix;

const contains = @import("helper.zig").contains;

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

        /// Clone the tensor. The function returns an error if the allocation of the new matrix fails.
        pub fn clone(self: Tensor(T)) !Tensor(T) {
            const other = try Tensor(T).init(self.shape, self.allocator);

            @memcpy(other.data, self.data);

            return other;
        }

        /// Fill the tensor with the specified value.
        pub fn fill(self: *Tensor(T), value: T) void {
            @memset(self.data, value);
        }

        pub fn matrix(self: Tensor(T), separator_index: usize) Matrix(T) {
            return Matrix(T){
                .data = self.data,
                .rows = mth.prod(usize, self.shape[0..separator_index]),
                .cols = mth.prod(usize, self.shape[separator_index.. ]),
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

        /// Fill the tensor with random values from a normal distribution with the specified mean and standard deviation.
        pub fn randn(self: *Tensor(T), mean: T, stdev: T, seed: u64) void {
            var rng = std.Random.DefaultPrng.init(seed); const random = rng.random();

            for (0..self.data.len) |i| {
                self.data[i] = mean + stdev * random.floatNorm(T);
            }
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

/// Contract two tensors using the TTGT scheme.
pub fn contract(comptime T: type, C: anytype, A: anytype, Ai: []const usize, B: anytype, Bi: []const usize, allocator: std.mem.Allocator) !void {
    const aisten = std.meta.fieldIndex(@TypeOf(A  ), "shape") != null;
    const bisten = std.meta.fieldIndex(@TypeOf(B  ), "shape") != null;
    const cisten = std.meta.fieldIndex(@TypeOf(C.*), "shape") != null;

    var a_shape = if (comptime aisten) A.shape else [2]usize{A.rows, A.cols};
    var b_shape = if (comptime bisten) B.shape else [2]usize{B.rows, B.cols};
    var c_shape = if (comptime cisten) C.shape else [2]usize{C.rows, C.cols};

    var a_stride = if (comptime aisten) A.stride else [2]usize{A.cols, 1};
    var b_stride = if (comptime bisten) B.stride else [2]usize{B.cols, 1};
    var c_stride = if (comptime cisten) C.stride else [2]usize{C.cols, 1};

    const A_TEN = if (comptime aisten) A   else Tensor(T){.data = A.data, .shape = &a_shape[0..].*, .stride = &a_stride[0..].*, .allocator = allocator};
    const B_TEN = if (comptime bisten) B   else Tensor(T){.data = B.data, .shape = &b_shape[0..].*, .stride = &b_stride[0..].*, .allocator = allocator};
    const C_TEN = if (comptime cisten) C.* else Tensor(T){.data = C.data, .shape = &c_shape[0..].*, .stride = &c_stride[0..].*, .allocator = allocator};

    var a_axes = std.ArrayList(usize).init(allocator); defer a_axes.deinit();
    var b_axes = std.ArrayList(usize).init(allocator); defer b_axes.deinit();

    for (0..A_TEN.shape.len) |i| if (!contains(usize, Ai, i)) {
        try a_axes.append(i);
    };

    for (Ai) |i| try a_axes.append(i);
    for (Bi) |i| try b_axes.append(i);

    for (0..B_TEN.shape.len) |i| if (!contains(usize, Bi, i)) {
        try b_axes.append(i);
    };

    const AT_shape = try allocator.alloc(usize, A_TEN.shape.len); defer allocator.free(AT_shape);
    const BT_shape = try allocator.alloc(usize, B_TEN.shape.len); defer allocator.free(BT_shape);

    for (0..A_TEN.shape.len) |i| AT_shape[i] = A_TEN.shape[a_axes.items[i]];
    for (0..B_TEN.shape.len) |i| BT_shape[i] = B_TEN.shape[b_axes.items[i]];

    var AT = try Tensor(T).init(AT_shape, allocator); defer AT.deinit();
    var BT = try Tensor(T).init(BT_shape, allocator); defer BT.deinit();

    try transpose(T, &AT, A_TEN, a_axes.items);
    try transpose(T, &BT, B_TEN, b_axes.items);

    const ATM = AT.matrix(a_axes.items.len - Ai.len); const BTM = BT.matrix(Bi.len);

    var CM = C_TEN.matrix(0); CM.rows = ATM.rows; CM.cols = BTM.cols;

    try cwp.Blas(T).dgemm(&CM, ATM, false, BTM, false);
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
pub fn transpose(comptime T: type, B: *Tensor(T), A: Tensor(T), axes: []const usize) !void {
    if (axes.len == 4) {
        for (0..A.shape[0]) |i| for (0..A.shape[1]) |j| for (0..A.shape[2]) |k| for (0..A.shape[3]) |l| {

            const source = [4]usize{i, j, k, l}; var destination: [4]usize = undefined;

            for (0..axes.len) |m| {
                destination[m] = source[axes[m]];
            }

            B.ptr(&destination).* = A.at(&source);
        };
    }

    else if (axes.len == 2) {
        for (0..A.shape[0]) |i| for (0..A.shape[1]) |j| {

            const source = [2]usize{i, j}; var destination: [2]usize = undefined;

            for (0..axes.len) |m| {
                destination[m] = source[axes[m]];
            }

            B.ptr(&destination).* = A.at(&source);
        };
    }

    else {
        return error.TensorTransposeNotImplemented;
    }
}
