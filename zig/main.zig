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

// fn Tensor(comptime T: type) type {
//     return struct {
//         data: []T,
//         shape: []i32,
//
//         pub fn debug(self: Tensor(T)) void {
//             return std.debug.print("{}", self.shape[1]);
//         }
//     };
// }

pub fn main() !void {
    const z = Complex(f64){ .re = 1.1, .im = 2.1 };

    // const A = Tensor(f64){ .data = [_]f64{ 1.1, 1.2, 1.3, 1.4 }, .shape = [_]i32{ 2, 2 } };

    z.debug();
    // A.debug();
}
