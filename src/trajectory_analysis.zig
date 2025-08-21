//! File with various trajectory analysis tools.

const std = @import("std"); const cwp = @import("cwrapper.zig");

const cnt = @import("constant.zig");
const hlp = @import("helper.zig"  );
const inp = @import("input.zig"   );
const out = @import("output.zig"  );

const Matrix       = @import("matrix.zig"       ).Matrix;
const StridedArray = @import("strided_array.zig").StridedArray;
const Vector       = @import("vector.zig"       ).Vector;

/// The main function of the trajectory module.
pub fn run(comptime T: type, opt: inp.TrajectoryAnalysisOptions(T), print: bool, allocator: std.mem.Allocator) !void {
    const parsed_trajectory = try parseTrajectory(T, opt.input, allocator);

    const atoms = parsed_trajectory.atoms; var coords = parsed_trajectory.coords; defer atoms.deinit(); defer coords.deinit();

    if (opt.remove_translation) removeTranslation(T, &coords);

    if (opt.remove_rotation) try removeRotation(T, &coords, allocator);

    if (opt.entropy != null) {

        var entropy: T = 0;

        if (std.mem.eql(u8, opt.entropy.?.algorithm, "schlitter")) entropy = try schlitterEntropy(T, coords, atoms, opt.entropy.?.temperature, allocator);

        if (print) try hlp.print(std.fs.File.stdout(), "\nSCHLITTER ENTROPY: {d:.14} a.u. = {d:.14} J/MOL\n", .{entropy, entropy * cnt.NA / cnt.J2AU});
    }

    if (opt.output) |path| try writeTrajectory(T, path, atoms, coords);
}

/// Parses the trajectory input file into a matrix, where each line is a frame and each column is a coordinate.
pub fn parseTrajectory(comptime T: type, path: []const u8, allocator: std.mem.Allocator) !struct {atoms: Vector(u32), coords: Matrix(T)} {
    const file = try std.fs.cwd().openFile(path, .{}); defer file.close();

    var buffer: [1024]u8 = undefined; var reader = file.reader(&buffer); var reader_interface = &reader.interface;

    var atoms  = std.ArrayList(u32){};
    var coords = std.ArrayList(T  ){};

    while (true) {

        const natom = try std.fmt.parseInt(usize, reader_interface.takeDelimiterExclusive('\n') catch break, 10);

        _ = try reader_interface.discardDelimiterInclusive('\n');

        for (0..natom) |_| {

            var it = std.mem.tokenizeAny(u8, try reader_interface.takeDelimiterExclusive('\n'), " ,");

            const symbol = it.next().?;

            if (atoms.items.len < natom) try atoms.append(allocator, cnt.SM2AN.get(symbol).?);

            const x = try std.fmt.parseFloat(T, it.next().?) * cnt.A2AU;
            const y = try std.fmt.parseFloat(T, it.next().?) * cnt.A2AU;
            const z = try std.fmt.parseFloat(T, it.next().?) * cnt.A2AU;

            try coords.append(allocator, x); try coords.append(allocator, y); try coords.append(allocator, z);
        }
    }

    const natom = atoms.items.len; const ntraj = coords.items.len / natom / 3;

    const coords_matrix = Matrix(T){.data = try coords.toOwnedSlice(allocator), .rows = ntraj, .cols = 3 * natom, .allocator = allocator};

    return .{.atoms = Vector(u32){.data = try atoms.toOwnedSlice(allocator), .rows = natom, .allocator = allocator}, .coords = coords_matrix};
}

/// Removes the rotation from the trajectory data by subtracting the mean position of each atom.
pub fn removeRotation(comptime T: type, coords: *Matrix(T), allocator: std.mem.Allocator) !void {
    var P = coords.row(0); try P.reshape(coords.cols / 3, 3); var QN = try P.clone();

    var H   = try Matrix(T).init(3, 3, allocator);
    var HU  = try Matrix(T).init(3, 3, allocator);
    var HS  = try Matrix(T).init(3, 3, allocator);
    var HVT = try Matrix(T).init(3, 3, allocator);
    var R   = try Matrix(T).init(3, 3, allocator);

    for (1..coords.rows) |i| {

        var Q = coords.row(i); try Q.reshape(coords.cols / 3, 3);

        try cwp.Blas(T).dgemm(&H, P, true, Q, false);

        try cwp.Lapack(T).dgesdd(&HU, &HS, &HVT, H);

        try cwp.Blas(T).dgemm(&R, HVT, true, HU, true);

        if (cwp.Eigen(T).det(R) < 0) {
            HVT.row(2).muls(-1); try cwp.Blas(T).dgemm(&R, HVT, true, HU, true);
        }

        try cwp.Blas(T).dgemm(&QN, Q, false, R, false); QN.memcpy(Q);
    }
}

/// Removes the translation from the trajectory data by subtracting the mean position of each atom.
pub fn removeTranslation(comptime T: type, coords: *Matrix(T)) void {
    for (0..coords.rows) |i| {

        const xarr = StridedArray(T){.data = coords.row(i).data, .len = coords.cols / 3, .stride = 3, .zero = 0};
        const yarr = StridedArray(T){.data = coords.row(i).data, .len = coords.cols / 3, .stride = 3, .zero = 1};
        const zarr = StridedArray(T){.data = coords.row(i).data, .len = coords.cols / 3, .stride = 3, .zero = 2};

        xarr.subs(xarr.mean()); yarr.subs(yarr.mean()); zarr.subs(zarr.mean());
    }
}

/// Function to calculate the Schlitter entropy of the trajectory data.
pub fn schlitterEntropy(comptime T: type, coords: Matrix(T), atoms: Vector(u32), temp: T, allocator: std.mem.Allocator) !T {
    var M = try Matrix(T).init(3 * atoms.rows, 3 * atoms.rows, allocator); defer M.deinit();
    var S = try Matrix(T).init(3 * atoms.rows, 3 * atoms.rows, allocator); defer S.deinit();

    for (0..3 * atoms.rows) |i| {
        M.ptr(i, i).* = cnt.AN2M[atoms.at(i / 3) - 1] * cnt.U2AU;
    }

    var sigma = try Matrix(T).init(3 * atoms.rows, 3 * atoms.rows, allocator); defer sigma.deinit();

    sigma.fill(0);

    for (0..sigma.rows) |i| for (0..sigma.cols) |j| {

        const xcol = coords.column(i);
        const ycol = coords.column(j);

        xcol.subs(xcol.mean());
        ycol.subs(ycol.mean());

        for (0..coords.rows) |k| {
            sigma.ptr(i, j).* += xcol.at(k) * ycol.at(k);
        }

        sigma.ptr(i, j).* /= hlp.asfloat(T, coords.rows - 1);
    };

    try cwp.Blas(T).dgemm(&S, M, false, sigma, false);

    S.muls(cnt.kB * cnt.J2AU * std.math.e * std.math.e * temp);

    for (0..S.rows) |i| S.ptr(i, i).* += 1;

    return 0.5 * cnt.kB * cnt.J2AU * std.math.log(T, std.math.e, cwp.Eigen(T).det(S));
}

/// Writes the trajectory data to the output file.
pub fn writeTrajectory(comptime T: type, path: []const u8, atoms: Vector(u32), coords: Matrix(T)) !void {
    const file = try std.fs.cwd().createFile(path, .{}); defer file.close();

    var buffer: [1024]u8 = undefined; var writer = file.writer(&buffer); var writer_interface = &writer.interface;

    for (0..coords.rows) |i| {

        try writer_interface.print("{d}\nGeometry #{d}\n", .{atoms.rows, i});

        for (0..atoms.rows) |j| {

            const x = coords.at(i, 3 * j + 0) / cnt.A2AU;
            const y = coords.at(i, 3 * j + 1) / cnt.A2AU;
            const z = coords.at(i, 3 * j + 2) / cnt.A2AU;

            try writer_interface.print("{s:2} {d:20.14} {d:20.14} {d:20.14}\n", .{try cnt.AN2SM(atoms.at(j)), x, y, z});
        }
    }

    try writer_interface.flush();
}
