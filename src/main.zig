const std = @import("std");
const debug = std.debug;
const math = std.math;
const ArrayList = std.ArrayList;

const rtw = @import("rtw.zig");

const Ray = rtw.ray.Ray;
const Vec3 = rtw.vec.Vec3;
const Point3 = rtw.vec.Point3;
const Color = rtw.vec.Color;
const rgb = rtw.vec.rgb;
const Random = rtw.rand.Random;
const randomPointInUnitSphere = rtw.rand.randomPointInUnitSphere;
const randomPointInUnitDisk = rtw.rand.randomPointInUnitDisk;
const randomUnitVector = rtw.rand.randomUnitVector;
const randomInt = rtw.rand.randomInt;
const randomReal01 = rtw.rand.randomReal01;
const randomReal = rtw.rand.randomReal;

const Material = rtw.material.Material;
const DiffuseMaterial = rtw.material.DiffuseMaterial;
const MetalMaterial = rtw.material.MetalMaterial;
const DielectricMaterial = rtw.material.DielectricMaterial;

const Texture = rtw.texture.Texture;
const SolidTexture = rtw.texture.SolidTexture;
const CheckerTexture = rtw.texture.CheckerTexture;

const Hittable = rtw.hittable.Hittable;
const Sphere = rtw.hittable.Sphere;
const MovingSphere = rtw.hittable.MovingSphere;
const HittableList = rtw.hittable.HittableList;
const BvhNode = rtw.hittable.BvhNode;

const Aabb = rtw.aabb.Aabb;
const HitRecord = rtw.hit_record.HitRecord;

fn makeSphere(center: Point3, radius: f64, material: *const Material) Hittable {
    return .{
        .sphere = .{
            .center = center,
            .radius = radius,
            .material = material,
        },
    };
}

const inf = math.inf(f64);
const pi = math.pi;

fn deg2rad(degree: f64) f64 {
    return degree * pi / 180.0;
}

const Camera = struct {
    origin: Point3,
    horizontal: Vec3,
    vertical: Vec3,
    lower_left_corner: Point3,
    u: Vec3,
    v: Vec3,
    w: Vec3,
    lens_radius: f64,
    time0: f64,
    time1: f64,

    fn init(
        lookFrom: Point3,
        lookAt: Point3,
        vup: Vec3,
        vfov: f64,
        aspect_ratio: f64,
        aperture: f64,
        focus_dist: f64,
        time0: f64,
        time1: f64,
    ) Camera {
        const theta = deg2rad(vfov);
        const h = @tan(theta / 2);
        const viewport_height = 2.0 * h;
        const viewport_width = aspect_ratio * viewport_height;

        const w = lookFrom.sub(lookAt).normalized();
        const u = vup.cross(w).normalized();
        const v = w.cross(u);

        const origin = lookFrom;
        const horizontal = u.mul(viewport_width * focus_dist);
        const vertical = v.mul(viewport_height * focus_dist);
        const lower_left_corner = origin.sub(horizontal.div(2.0)).sub(vertical.div(2.0)).sub(w.mul(focus_dist));

        return .{
            .origin = origin,
            .horizontal = horizontal,
            .vertical = vertical,
            .lower_left_corner = lower_left_corner,
            .u = u,
            .v = v,
            .w = w,
            .lens_radius = aperture / 2.0,
            .time0 = time0,
            .time1 = time1,
        };
    }

    fn getRay(camera: Camera, rng: Random, s: f64, t: f64) Ray {
        const rd = randomPointInUnitDisk(rng).mul(camera.lens_radius);
        const offset = camera.u.mul(rd.x).add(camera.v.mul(rd.y));
        const dir = camera.lower_left_corner.add(camera.horizontal.mul(s)).add(camera.vertical.mul(t)).sub(camera.origin).sub(offset);
        return .{
            .origin = camera.origin.add(offset),
            .dir = dir,
            .time = randomReal(rng, camera.time0, camera.time1),
        };
    }
};

fn rayColor(r: Ray, world: Hittable, rng: Random, depth: u32) Color {
    var rec: HitRecord = undefined;
    if (depth == 0) {
        // If we've exceeded the ray bounce limit, no more ligth is gathered.
        return rgb(0.0, 0.0, 0.0);
    }
    if (world.hit(r, 0.001, inf, &rec)) {
        var scattered: Ray = undefined;
        var attenuation: Color = undefined;
        if (rec.material.scatter(r, rec, &attenuation, &scattered, rng)) {
            return attenuation.mulV(rayColor(scattered, world, rng, depth - 1));
        } else {
            return rgb(0.0, 0.0, 0.0);
        }
    }
    const unit_dir = r.dir.normalized();
    const s = 0.5 * (unit_dir.y + 1.0);
    return rgb(1.0, 1.0, 1.0).mul(1.0 - s).add(rgb(0.5, 0.7, 1.0).mul(s));
}

fn writeColor(out: anytype, c: Color, samples_per_pixel: u32) !void {
    const scale = 1.0 / @intToFloat(f64, samples_per_pixel);
    try out.print("{} {} {}\n", .{
        @floatToInt(u8, 256.0 * math.clamp(@sqrt(c.x * scale), 0.0, 0.999)),
        @floatToInt(u8, 256.0 * math.clamp(@sqrt(c.y * scale), 0.0, 0.999)),
        @floatToInt(u8, 256.0 * math.clamp(@sqrt(c.z * scale), 0.0, 0.999)),
    });
}

fn generateTwoSpheres(rng: Random, allocator: anytype) !Hittable {
    _ = rng;
    var hittable_objects = ArrayList(Hittable).init(allocator);

    const checker = try Texture.makeChecker(allocator, rgb(0.2, 0.3, 0.1), rgb(0.9, 0.9, 0.9));
    var mat1 = try allocator.create(Material);
    var mat2 = try allocator.create(Material);

    mat1.* = .{ .diffuse = .{ .albedo = checker } };
    mat2.* = .{ .diffuse = .{ .albedo = checker } };

    try hittable_objects.append(makeSphere(.{ .x = 0, .y = -10, .z = 0 }, 10, mat1));
    try hittable_objects.append(makeSphere(.{ .x = 0, .y = 10, .z = 0 }, 10, mat2));

    return .{ .list = .{ .objects = hittable_objects } };
}

fn generateTwoPerlinSpheres(rng: Random, allocator: anytype) !Hittable {
    var hittable_objects = ArrayList(Hittable).init(allocator);

    const perlin = Texture.makeNoise(rng);
    var mat1 = try allocator.create(Material);
    var mat2 = try allocator.create(Material);

    mat1.* = .{ .diffuse = .{ .albedo = perlin } };
    mat2.* = .{ .diffuse = .{ .albedo = perlin } };

    try hittable_objects.append(makeSphere(.{ .x = 0, .y = -1000, .z = 0 }, 1000, mat1));
    try hittable_objects.append(makeSphere(.{ .x = 0, .y = 2, .z = 0 }, 2, mat2));

    return .{ .list = .{ .objects = hittable_objects } };
}

fn generateRandomScene(rng: Random, allocator: anytype) !Hittable {
    var hittable_objects = ArrayList(Hittable).init(allocator);

    var mat_ground = try allocator.create(Material);
    var mat1 = try allocator.create(Material);
    var mat2 = try allocator.create(Material);
    var mat3 = try allocator.create(Material);

    const checker = try Texture.makeChecker(allocator, rgb(0.2, 0.3, 0.1), rgb(0.9, 0.9, 0.9));

    mat_ground.* = .{ .diffuse = .{ .albedo = checker } };
    mat1.* = .{ .dielectric = .{ .ir = 1.5 } };
    mat2.* = .{ .diffuse = .{ .albedo = .{ .solid = .{ .color = rgb(0.4, 0.2, 0.1) } } } };
    mat3.* = .{ .metal = .{ .albedo = rgb(0.7, 0.6, 0.5), .fuzz = 0.0 } };

    try hittable_objects.append(makeSphere(.{ .x = 0, .y = -1000, .z = 0 }, 1000, mat_ground));
    try hittable_objects.append(makeSphere(.{ .x = 0, .y = 1, .z = 0 }, 1.0, mat1));
    try hittable_objects.append(makeSphere(.{ .x = -4, .y = 1, .z = 0 }, 1.0, mat2));
    try hittable_objects.append(makeSphere(.{ .x = 4, .y = 1, .z = 0 }, 1.0, mat3));

    var a: i32 = -3;
    while (a < 3) : (a += 1) {
        var b: i32 = -3;
        while (b < 3) : (b += 1) {
            const choose_mat = randomReal01(rng);
            const center = Point3{
                .x = @intToFloat(f64, a) + 0.9 * randomReal01(rng),
                .y = 0.2,
                .z = @intToFloat(f64, b) + 0.9 * randomReal01(rng),
            };

            if (center.sub(.{ .x = 4, .y = 0.2, .z = 0 }).norm() <= 0.9) {
                continue;
            }

            var mat_sphere = try allocator.create(Material);
            if (choose_mat < 0.8) {
                // diffuse
                const albedo = Color.random01(rng).mulV(Color.random01(rng));
                mat_sphere.* = .{ .diffuse = .{ .albedo = .{ .solid = .{ .color = albedo } } } };
                const center1 = center.add(.{ .x = 0, .y = randomReal(rng, 0, 0.5), .z = 0 });
                try hittable_objects.append(.{ .movingSphere = .{
                    .center0 = center,
                    .center1 = center1,
                    .time0 = 0,
                    .time1 = 1,
                    .radius = 0.2,
                    .material = mat_sphere,
                } });
            } else if (choose_mat < 0.95) {
                // metal
                const albedo = Color.random(rng, 0.5, 1);
                const fuzz = randomReal(rng, 0, 0.5);
                mat_sphere.* = .{ .metal = .{ .albedo = albedo, .fuzz = fuzz } };
                try hittable_objects.append(makeSphere(center, 0.2, mat_sphere));
            } else {
                // glass
                mat_sphere.* = .{ .dielectric = .{ .ir = 1.5 } };
                try hittable_objects.append(makeSphere(center, 0.2, mat_sphere));
            }
        }
    }

    return .{ .list = .{ .objects = hittable_objects } };
}

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const allocator = gpa.allocator();
    defer debug.assert(!gpa.deinit());

    var rng_ = std.rand.DefaultPrng.init(42);
    var rng = rng_.random();

    // Image
    const aspect_ratio = 3.0 / 2.0;
    const image_width = 600;
    const image_height = @floatToInt(comptime_int, @divTrunc(image_width, aspect_ratio));
    const samples_per_pixel = 50;
    const max_depth = 50;

    const scene = 3;

    // World
    var world: Hittable = undefined;
    var lookFrom: Point3 = undefined;
    var lookAt: Point3 = undefined;
    var vFov: f64 = 40;
    var aperture: f64 = 0;

    if (scene == 1) {
        world = try generateRandomScene(rng, allocator);
        lookFrom = .{ .x = 13, .y = 2, .z = 3 };
        lookAt = .{ .x = 0, .y = 0, .z = 0 };
        vFov = 20.0;
        aperture = 0.1;
    } else if (scene == 2) {
        world = try generateTwoSpheres(rng, allocator);
        lookFrom = .{ .x = 13, .y = 2, .z = 3 };
        lookAt = .{ .x = 0, .y = 0, .z = 0 };
        vFov = 20.0;
    } else if (scene == 3) {
        world = try generateTwoPerlinSpheres(rng, allocator);
        lookFrom = .{ .x = 13, .y = 2, .z = 3 };
        lookAt = .{ .x = 0, .y = 0, .z = 0 };
        vFov = 20.0;
    }
    defer world.deinit(allocator);

    // Camera
    const camera = Camera.init(
        lookFrom,
        lookAt,
        .{ .x = 0, .y = 1, .z = 0 },
        vFov,
        aspect_ratio,
        aperture,
        10.0,
        0,
        1,
    );

    // Render
    const stdout_file = std.io.getStdOut().writer();
    var bw = std.io.bufferedWriter(stdout_file);
    const stdout = bw.writer();

    try stdout.print("P3\n{} {}\n255\n", .{ image_width, image_height });

    var j: i32 = image_height - 1;
    while (j >= 0) : (j -= 1) {
        std.debug.print("\rScanlines remaining: {}     ", .{j});
        var i: i32 = 0;
        while (i < image_width) : (i += 1) {
            var s: u32 = 0;
            var pixelColor = rgb(0.0, 0.0, 0.0);
            while (s < samples_per_pixel) : (s += 1) {
                const u = (@intToFloat(f64, i) + randomReal01(rng)) / (image_width - 1);
                const v = (@intToFloat(f64, j) + randomReal01(rng)) / (image_height - 1);
                const r = camera.getRay(rng, u, v);
                pixelColor = pixelColor.add(rayColor(r, world, rng, max_depth));
            }
            try writeColor(stdout, pixelColor, samples_per_pixel);
        }
    }

    try bw.flush();
}
