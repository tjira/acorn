const std = @import("std");

fn Complex(comptime T: type) type {
    return struct {
        re: T,
        im: T,

        pub fn add(self: Complex(T), other: Complex(T)) Complex(T) {
            return Complex(T){
                .re = self.re + other.re,
                .im = self.im + other.im,
            };
        }
        pub fn debug(self: Complex(T)) void {
            return std.debug.print("({}, {})\n", .{ self.re, self.im });
        }
    };
}

fn Tensor(comptime T: type) type {
    return struct {
        data: std.ArrayList(T),
        shape: std.ArrayList(usize),

        pub fn debug(self: Tensor(T)) void {
            for (self.data.items) |item| std.debug.print("{}\n", .{item});
        }
        pub fn init() !Tensor(T) {
            return Tensor(T){ .data = std.ArrayList(T).init(std.heap.page_allocator), .shape = std.ArrayList(usize).init(std.heap.page_allocator) };
        }
        pub fn zero(shape: []const usize) !Tensor(T) {
            var tensor = try Tensor(T).init();
            try tensor.shape.appendSlice(shape);
            try tensor.data.ensureTotalCapacity(prod(usize, shape));
            for (prod(usize, shape)) |_| try tensor.data.append(0);
            return tensor;
        }
    };
}

fn prod(comptime T: type, array: []const T) T {
    var result: T = 1;

    for (array) |value| result *= value;

    return result;
}

pub fn main() !void {
    const A = try Tensor(f64).zero(&[_]usize{ 2, 4 });

    A.debug();
}
