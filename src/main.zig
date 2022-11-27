const std = @import("std");

const Vec3f = @Vector(3, f64);
const Point3 = Vec3f;
const Color = Vec3f;

fn vecNorm(v: Vec3f) f64 {
    return @sqrt(vecNormSquared(v));
}

fn vecNormSquared(v: Vec3f) f64 {
    return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}

fn vecDot(u: Vec3f, v: Vec3f) f64 {
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}

fn vecCross(u: Vec3f, v: Vec3f) Vec3f {
    return Vec3f(
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    );
}

fn vecNormalized(v: Vec3f) Vec3f {
    return v / vecNorm(v);
}

pub fn main() !void {
    // Image

    const image_width = 256;
    const image_height = 256;

    // Render

    const stdout_file = std.io.getStdOut().writer();
    var bw = std.io.bufferedWriter(stdout_file);
    const stdout = bw.writer();

    try stdout.print("P3\n{} {}\n255\n", .{ image_width, image_height });

    var j: i32 = image_height - 1;
    while (j >= 0) {
        std.debug.print("\rScanlines remaining: {}", .{j});
        var i: i32 = 0;
        while (i < image_width) {
            const r = @intToFloat(f64, i) / (image_width - 1);
            const g = @intToFloat(f64, j) / (image_height - 1);
            const b = 0.25;

            const ir = @floatToInt(u8, 255.999 * r);
            const ig = @floatToInt(u8, 255.999 * g);
            const ib = @floatToInt(u8, 255.999 * b);

            try stdout.print("{} {} {}\n", .{ ir, ig, ib });
            i += 1;
        }
        j -= 1;
    }

    try bw.flush();
}
