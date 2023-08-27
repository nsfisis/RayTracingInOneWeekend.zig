const std = @import("std");
const debug = std.debug;
const math = std.math;
const ArrayList = std.ArrayList;

const zigimg = @import("zigimg");
const Image = zigimg.Image;
const PixelFormat = zigimg.PixelFormat;

const Rc = @import("rc.zig").Rc;

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
        look_from: Point3,
        look_at: Point3,
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

        const w = look_from.sub(look_at).normalized();
        const u = vup.cross(w).normalized();
        const v = w.cross(u);

        const origin = look_from;
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

fn rayColor(r: Ray, background: Color, world: Hittable, rng: Random, depth: u32) Color {
    var rec: HitRecord = undefined;
    if (depth == 0) {
        // If we've exceeded the ray bounce limit, no more ligth is gathered.
        return rgb(0.0, 0.0, 0.0);
    }
    if (!world.hit(r, 0.001, inf, &rec)) {
        // If the ray hits nothing, return the background color.
        return background;
    }

    var scattered: Ray = undefined;
    var attenuation: Color = undefined;
    const emitted = rec.material.emitted(rec.u, rec.v, rec.p);
    if (rec.material.scatter(r, rec, &attenuation, &scattered, rng)) {
        return emitted.add(attenuation.mulV(rayColor(scattered, background, world, rng, depth - 1)));
    } else {
        return emitted;
    }
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

    const perlin = try Texture.makeNoise(allocator, 4.0, rng);
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
                .x = @as(f64, @floatFromInt(a)) + 0.9 * randomReal01(rng),
                .y = 0.2,
                .z = @as(f64, @floatFromInt(b)) + 0.9 * randomReal01(rng),
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

fn generateEarthScene(allocator: anytype) !Hittable {
    var hittable_objects = ArrayList(Hittable).init(allocator);

    const earth_texture = try Texture.makeImage(allocator, "assets/sekaichizu.png");
    var earth_surface = try allocator.create(Material);

    earth_surface.* = .{ .diffuse = .{ .albedo = earth_texture } };

    try hittable_objects.append(makeSphere(.{ .x = 0, .y = 0, .z = 0 }, 2, earth_surface));

    return .{ .list = .{ .objects = hittable_objects } };
}

fn generateSimpleLightScene(rng: Random, allocator: anytype) !Hittable {
    var hittable_objects = ArrayList(Hittable).init(allocator);

    const perlin = try Texture.makeNoise(allocator, 4.0, rng);
    var mat1 = try allocator.create(Material);
    var mat2 = try allocator.create(Material);

    mat1.* = .{ .diffuse = .{ .albedo = perlin } };
    mat2.* = .{ .diffuse = .{ .albedo = perlin } };

    try hittable_objects.append(makeSphere(.{ .x = 0, .y = -1000, .z = 0 }, 1000, mat1));
    try hittable_objects.append(makeSphere(.{ .x = 0, .y = 2, .z = 0 }, 2, mat2));

    const light = Texture.makeSolid(rgb(4, 4, 4));
    var mat3 = try Rc(Material).init(allocator);

    mat3.get_mut().* = .{ .diffuse_light = .{ .emit = light } };

    try hittable_objects.append(.{ .xyRect = .{ .x0 = 3.0, .x1 = 5.0, .y0 = 1.0, .y1 = 3.0, .k = -2.0, .material = mat3 } });

    return .{ .list = .{ .objects = hittable_objects } };
}

fn generateCornellBox(allocator: anytype) !Hittable {
    var hittable_objects = ArrayList(Hittable).init(allocator);

    var red = try Rc(Material).init(allocator);
    var white = try Rc(Material).init(allocator);
    var green = try Rc(Material).init(allocator);
    var light = try Rc(Material).init(allocator);

    red.get_mut().* = .{ .diffuse = .{ .albedo = Texture.makeSolid(rgb(0.65, 0.05, 0.05)) } };
    white.get_mut().* = .{ .diffuse = .{ .albedo = Texture.makeSolid(rgb(0.73, 0.73, 0.73)) } };
    green.get_mut().* = .{ .diffuse = .{ .albedo = Texture.makeSolid(rgb(0.12, 0.45, 0.15)) } };
    light.get_mut().* = .{ .diffuse_light = .{ .emit = Texture.makeSolid(rgb(15, 15, 15)) } };

    try hittable_objects.append(.{ .yzRect = .{ .y0 = 0, .y1 = 555, .z0 = 0, .z1 = 555, .k = 555, .material = green } });
    try hittable_objects.append(.{ .yzRect = .{ .y0 = 0, .y1 = 555, .z0 = 0, .z1 = 555, .k = 0, .material = red } });
    try hittable_objects.append(.{ .xzRect = .{ .x0 = 213, .x1 = 343, .z0 = 227, .z1 = 332, .k = 554, .material = light } });
    try hittable_objects.append(.{ .xzRect = .{ .x0 = 0, .x1 = 555, .z0 = 0, .z1 = 555, .k = 0, .material = white } });
    try hittable_objects.append(.{ .xzRect = .{ .x0 = 0, .x1 = 555, .z0 = 0, .z1 = 555, .k = 555, .material = white.clone() } });
    try hittable_objects.append(.{ .xyRect = .{ .x0 = 0, .x1 = 555, .y0 = 0, .y1 = 555, .k = 555, .material = white.clone() } });

    var box1 = try Rc(Hittable).init(allocator);
    var box1r = try Rc(Hittable).init(allocator);
    var box2 = try Rc(Hittable).init(allocator);
    var box2r = try Rc(Hittable).init(allocator);

    box1.get_mut().* = try Hittable.makeBox(.{ .x = 0, .y = 0, .z = 0 }, .{ .x = 165, .y = 330, .z = 165 }, white.clone(), allocator);
    box1r.get_mut().* = Hittable.rotateY(box1, deg2rad(15));
    try hittable_objects.append(Hittable.translate(box1r, .{ .x = 265, .y = 0, .z = 295 }));

    box2.get_mut().* = try Hittable.makeBox(.{ .x = 0, .y = 0, .z = 0 }, .{ .x = 165, .y = 165, .z = 165 }, white.clone(), allocator);
    box2r.get_mut().* = Hittable.rotateY(box2, deg2rad(-18));
    try hittable_objects.append(Hittable.translate(box2r, .{ .x = 130, .y = 0, .z = 65 }));

    return .{ .list = .{ .objects = hittable_objects } };
}

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const allocator = gpa.allocator();
    defer debug.assert(gpa.deinit() == .ok);

    var rng_ = std.rand.DefaultPrng.init(42);
    var rng = rng_.random();

    // Image
    var aspect_ratio: f64 = 3.0 / 2.0;
    var image_width: u32 = 600;
    var image_height: u32 = @as(u32, @intFromFloat(@divTrunc(@as(f64, @floatFromInt(image_width)), aspect_ratio)));
    const max_depth = 50;
    var samples_per_pixel: u32 = 50;

    const scene = 6;

    // World
    var world: Hittable = undefined;
    var look_from: Point3 = undefined;
    var look_at: Point3 = undefined;
    var vfov: f64 = 40;
    var aperture: f64 = 0;
    var background: Color = undefined;

    if (scene == 1) {
        world = try generateRandomScene(rng, allocator);
        background = rgb(0.70, 0.80, 1.00);
        look_from = .{ .x = 13, .y = 2, .z = 3 };
        look_at = .{ .x = 0, .y = 0, .z = 0 };
        vfov = 20.0;
        aperture = 0.1;
    } else if (scene == 2) {
        world = try generateTwoSpheres(rng, allocator);
        background = rgb(0.70, 0.80, 1.00);
        look_from = .{ .x = 13, .y = 2, .z = 3 };
        look_at = .{ .x = 0, .y = 0, .z = 0 };
        vfov = 20.0;
    } else if (scene == 3) {
        world = try generateTwoPerlinSpheres(rng, allocator);
        background = rgb(0.70, 0.80, 1.00);
        look_from = .{ .x = 13, .y = 2, .z = 3 };
        look_at = .{ .x = 0, .y = 0, .z = 0 };
        vfov = 20.0;
    } else if (scene == 4) {
        world = try generateEarthScene(allocator);
        background = rgb(0.70, 0.80, 1.00);
        look_from = .{ .x = 13, .y = 2, .z = 3 };
        look_at = .{ .x = 0, .y = 0, .z = 0 };
        vfov = 20.0;
    } else if (scene == 5) {
        world = try generateSimpleLightScene(rng, allocator);
        background = rgb(0, 0, 0);
        look_from = .{ .x = 26, .y = 3, .z = 6 };
        look_at = .{ .x = 0, .y = 2, .z = 0 };
        vfov = 20.0;
        samples_per_pixel = 400;
    } else if (scene == 6) {
        world = try generateCornellBox(allocator);
        background = rgb(0, 0, 0);
        look_from = .{ .x = 278, .y = 278, .z = -800 };
        look_at = .{ .x = 278, .y = 278, .z = 0 };
        vfov = 40.0;
        aspect_ratio = 1.0;
        image_width = 600;
        image_height = 600;
        samples_per_pixel = 200;
    }
    defer world.deinit();

    // Camera
    const camera = Camera.init(
        look_from,
        look_at,
        .{ .x = 0, .y = 1, .z = 0 },
        vfov,
        aspect_ratio,
        aperture,
        10.0,
        0,
        1,
    );

    // Render
    var image = try Image.create(allocator, image_width, image_height, PixelFormat.rgb24);
    defer image.deinit();

    var j: u32 = 0;
    while (j < image_height) : (j += 1) {
        std.debug.print("\rScanlines remaining: {}     ", .{image_height - j - 1});
        var i: u32 = 0;
        while (i < image_width) : (i += 1) {
            var s: u32 = 0;
            var pixelColor = rgb(0.0, 0.0, 0.0);
            while (s < samples_per_pixel) : (s += 1) {
                const u = (@as(f64, @floatFromInt(i)) + randomReal01(rng)) / (@as(f64, @floatFromInt(image_width)) - 1.0);
                const v = (@as(f64, @floatFromInt(j)) + randomReal01(rng)) / (@as(f64, @floatFromInt(image_height)) - 1.0);
                const r = camera.getRay(rng, u, v);
                pixelColor = pixelColor.add(rayColor(r, background, world, rng, max_depth));
            }
            const scale = 1.0 / @as(f64, @floatFromInt(samples_per_pixel));
            image.pixels.rgb24[i + (image_height - j - 1) * image_width] = .{
                .r = @intFromFloat(256.0 * math.clamp(@sqrt(pixelColor.x * scale), 0.0, 0.999)),
                .g = @intFromFloat(256.0 * math.clamp(@sqrt(pixelColor.y * scale), 0.0, 0.999)),
                .b = @intFromFloat(256.0 * math.clamp(@sqrt(pixelColor.z * scale), 0.0, 0.999)),
            };
        }
    }

    // Output
    try image.writeToFilePath("out.png", .{ .png = .{} });
}
