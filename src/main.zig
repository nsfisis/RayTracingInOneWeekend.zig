const std = @import("std");

const Vec3 = struct {
    x: f64,
    y: f64,
    z: f64,

    pub fn norm(v: Vec3) f64 {
        return @sqrt(v.normSquared());
    }

    pub fn normSquared(v: Vec3) f64 {
        return v.x * v.x + v.y * v.y + v.z * v.z;
    }

    pub fn dot(u: Vec3, v: Vec3) f64 {
        return u.x * v.x + u.y * v.y + u.z * v.z;
    }

    pub fn cross(u: Vec3, v: Vec3) Vec3 {
        return Vec3{
            .x = u.y * v.z - u.z * v.y,
            .y = u.z * v.x - u.x * v.z,
            .z = u.x * v.y - u.y * v.x,
        };
    }

    pub fn normalized(v: Vec3) Vec3 {
        return v.div(v.norm());
    }

    pub fn div(v: Vec3, t: f64) Vec3 {
        return Vec3{
            .x = v.x / t,
            .y = v.y / t,
            .z = v.z / t,
        };
    }
};

const Point3 = Vec3;
const Color = Vec3;

fn writeColor(out: anytype, c: Color) !void {
    try out.print("{} {} {}\n", .{
        @floatToInt(u8, 255.999 * c.x),
        @floatToInt(u8, 255.999 * c.y),
        @floatToInt(u8, 255.999 * c.z),
    });
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
    while (j >= 0) : (j -= 1) {
        std.debug.print("\rScanlines remaining: {}", .{j});
        var i: i32 = 0;
        while (i < image_width) : (i += 1) {
            try writeColor(stdout, Color{
                .x = @intToFloat(f64, i) / (image_width - 1),
                .y = @intToFloat(f64, j) / (image_height - 1),
                .z = 0.25,
            });
        }
    }

    try bw.flush();
}
