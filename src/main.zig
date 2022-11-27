const std = @import("std");
const debug = std.debug;

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
        const n = v.norm();
        if (n == 0.0) {
            return v;
        } else {
            return v.div(n);
        }
    }

    pub fn add(u: Vec3, v: Vec3) Vec3 {
        return Vec3{
            .x = u.x + v.x,
            .y = u.y + v.y,
            .z = u.z + v.z,
        };
    }

    pub fn sub(u: Vec3, v: Vec3) Vec3 {
        return Vec3{
            .x = u.x - v.x,
            .y = u.y - v.y,
            .z = u.z - v.z,
        };
    }

    pub fn mul(v: Vec3, t: f64) Vec3 {
        return Vec3{
            .x = v.x * t,
            .y = v.y * t,
            .z = v.z * t,
        };
    }

    pub fn div(v: Vec3, t: f64) Vec3 {
        return Vec3{
            .x = v.x / t,
            .y = v.y / t,
            .z = v.z / t,
        };
    }

    pub fn pp(v: Vec3) void {
        debug.print("{} {} {}\n", .{ v.x, v.y, v.z });
    }
};

const Point3 = Vec3;
const Color = Vec3;

const Ray = struct {
    origin: Vec3,
    dir: Vec3,

    pub fn at(r: Ray, t: f64) Point3 {
        return r.origin.add(r.dir.mul(t));
    }
};

fn hitSphere(center: Point3, radius: f64, r: Ray) f64 {
    const oc = r.origin.sub(center);
    const a = Vec3.dot(r.dir, r.dir);
    const b = 2.0 * Vec3.dot(oc, r.dir);
    const c = Vec3.dot(oc, oc) - radius * radius;
    const discriminant = b * b - 4.0 * a * c;
    if (discriminant < 0.0) {
        // r does not intersect the sphere.
        return -1.0;
    } else {
        // root of the equation.
        return (-b - @sqrt(discriminant)) / (2.0 * a);
    }
}

fn rayColor(r: Ray) Color {
    const sphereOrigin = Point3{ .x = 0.0, .y = 0.0, .z = -1.0 };
    const sphereRadius = 0.5;

    const t = hitSphere(sphereOrigin, sphereRadius, r);
    if (t > 0.0) {
        // r hints the sphere.
        const sphereNormal = r.at(t).sub(sphereOrigin).normalized();
        return (Color{ .x = sphereNormal.x + 1.0, .y = sphereNormal.y + 1.0, .z = sphereNormal.z + 1.0 }).mul(0.5);
    }
    const unit_dir = r.dir.normalized();
    const s = 0.5 * (unit_dir.y + 1.0);
    return (Color{ .x = 1.0, .y = 1.0, .z = 1.0 }).mul(1.0 - s).add((Color{ .x = 0.5, .y = 0.7, .z = 1.0 }).mul(s));
}

fn writeColor(out: anytype, c: Color) !void {
    try out.print("{} {} {}\n", .{
        @floatToInt(u8, 255.999 * c.x),
        @floatToInt(u8, 255.999 * c.y),
        @floatToInt(u8, 255.999 * c.z),
    });
}

pub fn main() !void {
    // Image

    const aspect_ratio = 16.0 / 9.0;
    const image_width = 400;
    const image_height = @floatToInt(comptime_int, @divTrunc(image_width, aspect_ratio));

    // Camera
    const viewport_height = 2.0;
    const viewport_width = aspect_ratio * viewport_height;
    const focal_length = 1.0;

    const origin = Point3{ .x = 0.0, .y = 0.0, .z = 0.0 };
    const horizontal = Vec3{ .x = viewport_width, .y = 0.0, .z = 0.0 };
    const vertial = Vec3{ .x = 0.0, .y = viewport_height, .z = 0.0 };
    const lower_left_corner = origin.sub(horizontal.div(2.0)).sub(vertial.div(2.0)).sub(Vec3{ .x = 0.0, .y = 0.0, .z = focal_length });

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
            const u = @intToFloat(f64, i) / (image_width - 1);
            const v = @intToFloat(f64, j) / (image_height - 1);
            const dir = lower_left_corner.add(horizontal.mul(u)).add(vertial.mul(v)).sub(origin);
            const r = Ray{ .origin = origin, .dir = dir };
            const pixelColor = rayColor(r);
            try writeColor(stdout, pixelColor);
        }
    }

    try bw.flush();
}
