const std = @import("std");
const debug = std.debug;
const math = std.math;
const ArrayList = std.ArrayList;

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

    pub fn mulV(u: Vec3, v: Vec3) Vec3 {
        return Vec3{
            .x = u.x * v.x,
            .y = u.y * v.y,
            .z = u.z * v.z,
        };
    }

    pub fn div(v: Vec3, t: f64) Vec3 {
        return Vec3{
            .x = v.x / t,
            .y = v.y / t,
            .z = v.z / t,
        };
    }

    pub fn random01(rand: std.rand.Random) Vec3 {
        return Vec3{
            .x = randomReal01(rand),
            .y = randomReal01(rand),
            .z = randomReal01(rand),
        };
    }

    pub fn random(rand: std.rand.Random, min: f64, max: f64) Vec3 {
        return Vec3{
            .x = randomReal(rand, min, max),
            .y = randomReal(rand, min, max),
            .z = randomReal(rand, min, max),
        };
    }

    // for debugging
    pub fn pp(v: Vec3) void {
        debug.print("{} {} {}\n", .{ v.x, v.y, v.z });
    }

    fn near_zero(v: Vec3) bool {
        const epsilon = 1e-8;
        return @fabs(v.x) < epsilon and @fabs(v.y) < epsilon and @fabs(v.z) < epsilon;
    }
};

fn randomPointInUnitSphere(rand: std.rand.Random) Vec3 {
    while (true) {
        const p = Vec3.random(rand, -1.0, 1.0);
        if (p.norm() >= 1) continue;
        return p;
    }
}

fn randomPointInUnitDisk(rand: std.rand.Random) Vec3 {
    while (true) {
        const p = Vec3{ .x = randomReal(rand, -1.0, 1.0), .y = randomReal(rand, -1.0, 1.0), .z = 0.0 };
        if (p.norm() >= 1) continue;
        return p;
    }
}

fn randomUnitVector(rand: std.rand.Random) Vec3 {
    return randomPointInUnitSphere(rand).normalized();
}

fn reflect(v: Vec3, n: Vec3) Vec3 {
    return v.sub(n.mul(2 * v.dot(n)));
}

fn refract(uv: Vec3, n: Vec3, etai_over_etat: f64) Vec3 {
    const cos_theta = @min(uv.mul(-1.0).dot(n), 1.0);
    const r_out_perpendicular = uv.add(n.mul(cos_theta)).mul(etai_over_etat);
    const r_out_parallel = n.mul(-@sqrt(@fabs(1.0 - r_out_perpendicular.normSquared())));
    return r_out_perpendicular.add(r_out_parallel);
}

const Point3 = Vec3;
const Color = Vec3;

fn rgb(r: f64, g: f64, b: f64) Color {
    return Color{ .x = r, .y = g, .z = b };
}

const Ray = struct {
    origin: Vec3,
    dir: Vec3,

    pub fn at(r: Ray, t: f64) Point3 {
        return r.origin.add(r.dir.mul(t));
    }
};

const MaterialTag = enum {
    diffuse,
    metal,
    dielectric,
};

const Material = union(MaterialTag) {
    diffuse: DiffuseMaterial,
    metal: MetalMaterial,
    dielectric: DielectricMaterial,

    fn scatter(mat: Material, r_in: Ray, record: HitRecord, attenuation: *Color, scattered: *Ray, rand: std.rand.Random) bool {
        return switch (mat) {
            MaterialTag.diffuse => |diffuse_mat| diffuse_mat.scatter(r_in, record, attenuation, scattered, rand),
            MaterialTag.metal => |metal_mat| metal_mat.scatter(r_in, record, attenuation, scattered, rand),
            MaterialTag.dielectric => |dielectric_mat| dielectric_mat.scatter(r_in, record, attenuation, scattered, rand),
        };
    }
};

const DiffuseMaterial = struct {
    albedo: Color,

    fn scatter(mat: DiffuseMaterial, r_in: Ray, record: HitRecord, attenuation: *Color, scattered: *Ray, rand: std.rand.Random) bool {
        _ = r_in;
        var scatter_direction = record.normal.add(randomUnitVector(rand));
        if (scatter_direction.near_zero()) {
            scatter_direction = record.normal;
        }
        scattered.* = Ray{ .origin = record.p, .dir = scatter_direction };
        attenuation.* = mat.albedo;
        return true;
    }
};

const MetalMaterial = struct {
    albedo: Color,
    fuzz: f64,

    fn scatter(mat: MetalMaterial, r_in: Ray, record: HitRecord, attenuation: *Color, scattered: *Ray, rand: std.rand.Random) bool {
        debug.assert(mat.fuzz <= 1.0);
        const reflected = reflect(r_in.dir.normalized(), record.normal);
        scattered.* = Ray{ .origin = record.p, .dir = reflected.add(randomPointInUnitSphere(rand).mul(mat.fuzz)) };
        attenuation.* = mat.albedo;
        return reflected.dot(record.normal) > 0.0;
    }
};

const DielectricMaterial = struct {
    // index of refraction.
    ir: f64,

    fn scatter(mat: DielectricMaterial, r_in: Ray, record: HitRecord, attenuation: *Color, scattered: *Ray, rand: std.rand.Random) bool {
        const refraction_ratio = if (record.front_face) 1.0 / mat.ir else mat.ir;
        const unit_dir = r_in.dir.normalized();

        const cos_theta = @min(unit_dir.mul(-1.0).dot(record.normal), 1.0);
        const sin_theta = @sqrt(1.0 - cos_theta * cos_theta);
        // sin(theta') = refraction_ratio * sin(theta) <= 1
        const can_refract = refraction_ratio * sin_theta <= 1.0;

        const dir = if (can_refract and reflectance(cos_theta, refraction_ratio) < randomReal01(rand)) refract(unit_dir, record.normal, refraction_ratio) else reflect(unit_dir, record.normal);
        scattered.* = Ray{ .origin = record.p, .dir = dir };
        attenuation.* = rgb(1.0, 1.0, 1.0);
        return true;
    }

    fn reflectance(cos: f64, refraction_idx: f64) f64 {
        const r0 = (1.0 - refraction_idx) / (1.0 + refraction_idx);
        const r1 = r0 * r0;
        return r1 + (1.0 - r1) * math.pow(f64, 1.0 - cos, 5.0);
    }
};

const HitRecord = struct {
    // The point where the ray and the hittable hits.
    p: Point3,
    // The normal of the hittable at p.
    normal: Vec3,
    // The material at p.
    material: *const Material,
    // p = ray.at(t)
    t: f64,
    // True if the ray hits the hittable from the front face, i.e., outside of it.
    front_face: bool,
};

const HittableTag = enum {
    sphere,
    list,
};

const Hittable = union(HittableTag) {
    sphere: Sphere,
    list: HittableList,

    fn hit(h: Hittable, r: Ray, t_min: f64, t_max: f64, record: *HitRecord) bool {
        return switch (h) {
            HittableTag.sphere => |sphere| sphere.hit(r, t_min, t_max, record),
            HittableTag.list => |list| list.hit(r, t_min, t_max, record),
        };
    }
};

const Sphere = struct {
    center: Point3,
    radius: f64,
    material: *const Material,

    fn hit(sphere: Sphere, r: Ray, t_min: f64, t_max: f64, record: *HitRecord) bool {
        const oc = r.origin.sub(sphere.center);
        const a = r.dir.normSquared();
        const half_b = Vec3.dot(oc, r.dir);
        const c = oc.normSquared() - sphere.radius * sphere.radius;

        const discriminant = half_b * half_b - a * c;
        if (discriminant < 0.0) {
            // r does not intersect the sphere.
            return false;
        }
        const sqrtd = @sqrt(discriminant);

        // Find the nearest root that lies in the acceptable range.
        var root = (-half_b - sqrtd) / a;
        if (root < t_min or t_max < root) {
            root = (-half_b + sqrtd) / a;
            if (root < t_min or t_max < root) {
                // out of range
                return false;
            }
        }

        record.t = root;
        record.p = r.at(root);
        const outward_normal = (record.p.sub(sphere.center)).div(sphere.radius);
        record.front_face = Vec3.dot(outward_normal, r.dir) < 0.0;
        if (record.front_face) {
            record.normal = outward_normal;
        } else {
            record.normal = outward_normal.mul(-1.0);
        }
        record.material = sphere.material;

        return true;
    }
};

const HittableList = struct {
    objects: ArrayList(*const Hittable),

    fn hit(list: HittableList, r: Ray, t_min: f64, t_max: f64, record: *HitRecord) bool {
        var hit_anything = false;
        var closest_so_far = t_max;

        for (list.objects.items) |object| {
            var rec: HitRecord = undefined;
            if (object.hit(r, t_min, closest_so_far, &rec)) {
                hit_anything = true;
                closest_so_far = rec.t;
                record.* = rec;
            }
        }
        return hit_anything;
    }
};

const inf = math.inf(f64);
const pi = math.pi;

fn deg2rad(degree: f64) f64 {
    return degree * pi / 180.0;
}

// [0, 1)
fn randomReal01(rand: std.rand.Random) f64 {
    return rand.float(f64);
}

// [min, max)
fn randomReal(rand: std.rand.Random, min: f64, max: f64) f64 {
    return min + randomReal01(rand) * (max - min);
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

    fn init(lookFrom: Point3, lookAt: Point3, vup: Vec3, vfov: f64, aspect_ratio: f64, aperture: f64, focus_dist: f64) Camera {
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

        return Camera{
            .origin = origin,
            .horizontal = horizontal,
            .vertical = vertical,
            .lower_left_corner = lower_left_corner,
            .u = u,
            .v = v,
            .w = w,
            .lens_radius = aperture / 2.0,
        };
    }

    fn getRay(camera: Camera, rand: std.rand.Random, s: f64, t: f64) Ray {
        const rd = randomPointInUnitDisk(rand).mul(camera.lens_radius);
        const offset = camera.u.mul(rd.x).add(camera.v.mul(rd.y));
        const dir = camera.lower_left_corner.add(camera.horizontal.mul(s)).add(camera.vertical.mul(t)).sub(camera.origin).sub(offset);
        return Ray{
            .origin = camera.origin.add(offset),
            .dir = dir,
        };
    }
};

fn rayColor(r: Ray, world: Hittable, rand: std.rand.Random, depth: u32) Color {
    var rec: HitRecord = undefined;
    if (depth == 0) {
        // If we've exceeded the ray bounce limit, no more ligth is gathered.
        return rgb(0.0, 0.0, 0.0);
    }
    if (world.hit(r, 0.001, inf, &rec)) {
        var scattered: Ray = undefined;
        var attenuation: Color = undefined;
        if (rec.material.scatter(r, rec, &attenuation, &scattered, rand)) {
            return attenuation.mulV(rayColor(scattered, world, rand, depth - 1));
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

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const allocator = gpa.allocator();
    defer debug.assert(!gpa.deinit());

    var rng = std.rand.DefaultPrng.init(42);
    var rand = rng.random();

    // Image
    const aspect_ratio = 16.0 / 9.0;
    const image_width = 400;
    const image_height = @floatToInt(comptime_int, @divTrunc(image_width, aspect_ratio));
    const samples_per_pixel = 100;
    const max_depth = 50;

    const material_ground = Material{ .diffuse = DiffuseMaterial{ .albedo = rgb(0.8, 0.8, 0.0) } };
    const material_center = Material{ .diffuse = DiffuseMaterial{ .albedo = rgb(0.1, 0.2, 0.5) } };
    const material_left = Material{ .dielectric = DielectricMaterial{ .ir = 1.5 } };
    const material_right = Material{ .metal = MetalMaterial{ .albedo = rgb(0.8, 0.6, 0.2), .fuzz = 0.0 } };
    const sphere1 = Hittable{ .sphere = Sphere{ .center = Point3{ .x = 0.0, .y = -100.5, .z = -1.0 }, .radius = 100.0, .material = &material_ground } };
    const sphere2 = Hittable{ .sphere = Sphere{ .center = Point3{ .x = 0.0, .y = 0.0, .z = -1.0 }, .radius = 0.5, .material = &material_center } };
    const sphere3 = Hittable{ .sphere = Sphere{ .center = Point3{ .x = -1.0, .y = 0.0, .z = -1.0 }, .radius = 0.5, .material = &material_left } };
    const sphere4 = Hittable{ .sphere = Sphere{ .center = Point3{ .x = -1.0, .y = 0.0, .z = -1.0 }, .radius = -0.45, .material = &material_left } };
    const sphere5 = Hittable{ .sphere = Sphere{ .center = Point3{ .x = 1.0, .y = 0.0, .z = -1.0 }, .radius = 0.5, .material = &material_right } };
    var hittable_objects = ArrayList(*const Hittable).init(allocator);
    try hittable_objects.append(&sphere1);
    try hittable_objects.append(&sphere2);
    try hittable_objects.append(&sphere3);
    try hittable_objects.append(&sphere4);
    try hittable_objects.append(&sphere5);
    const world = Hittable{ .list = HittableList{ .objects = hittable_objects } };
    defer hittable_objects.deinit();

    // Camera
    const lookFrom = Point3{ .x = 3.0, .y = 3.0, .z = 2.0 };
    const lookAt = Point3{ .x = 0.0, .y = 0.0, .z = -1.0 };
    const camera = Camera.init(
        lookFrom,
        lookAt,
        Vec3{ .x = 0.0, .y = 1.0, .z = 0.0 },
        20.0,
        aspect_ratio,
        2.0,
        lookFrom.sub(lookAt).norm(),
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
                const u = (@intToFloat(f64, i) + randomReal01(rand)) / (image_width - 1);
                const v = (@intToFloat(f64, j) + randomReal01(rand)) / (image_height - 1);
                const r = camera.getRay(rand, u, v);
                pixelColor = pixelColor.add(rayColor(r, world, rand, max_depth));
            }
            try writeColor(stdout, pixelColor, samples_per_pixel);
        }
    }

    try bw.flush();
}
