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
    return .{ .x = r, .y = g, .z = b };
}

const Ray = struct {
    origin: Vec3,
    dir: Vec3,
    time: f64,

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
        var scatter_direction = record.normal.add(randomUnitVector(rand));
        if (scatter_direction.near_zero()) {
            scatter_direction = record.normal;
        }
        scattered.* = .{ .origin = record.p, .dir = scatter_direction, .time = r_in.time };
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
        scattered.* = .{ .origin = record.p, .dir = reflected.add(randomPointInUnitSphere(rand).mul(mat.fuzz)), .time = r_in.time };
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
        scattered.* = .{ .origin = record.p, .dir = dir, .time = r_in.time };
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
    // The coordinate of the surface where the ray intersects.
    u: f64,
    v: f64,
    // True if the ray hits the hittable from the front face, i.e., outside of it.
    front_face: bool,
};

const HittableTag = enum {
    sphere,
    movingSphere,
    list,
    bvhNode,
};

const Hittable = union(HittableTag) {
    sphere: Sphere,
    movingSphere: MovingSphere,
    list: HittableList,
    bvhNode: BvhNode,

    fn hit(h: Hittable, r: Ray, tMin: f64, tMax: f64, record: *HitRecord) bool {
        return switch (h) {
            HittableTag.sphere => |sphere| sphere.hit(r, tMin, tMax, record),
            HittableTag.movingSphere => |movingSphere| movingSphere.hit(r, tMin, tMax, record),
            HittableTag.list => |list| list.hit(r, tMin, tMax, record),
            HittableTag.bvhNode => |node| node.hit(r, tMin, tMax, record),
        };
    }

    fn boudingBox(h: Hittable, time0: f64, time1: f64, outputBox: *Aabb) bool {
        return switch (h.*) {
            HittableTag.sphere => |sphere| sphere.boudingBox(time0, time1, outputBox),
            HittableTag.movingSphere => |movingSphere| movingSphere.boudingBox(time0, time1, outputBox),
            HittableTag.list => |list| list.boudingBox(time0, time1, outputBox),
            HittableTag.bvhNode => |node| node.boudingBox(time0, time1, outputBox),
        };
    }

    fn deinit(h: *const Hittable, allocator: anytype) void {
        return switch (h.*) {
            HittableTag.sphere => |sphere| sphere.deinit(allocator),
            HittableTag.movingSphere => |movingSphere| movingSphere.deinit(allocator),
            HittableTag.list => |list| list.deinit(allocator),
            HittableTag.bvhNode => |node| node.deinit(allocator),
        };
    }
};

const Sphere = struct {
    center: Point3,
    radius: f64,
    material: *const Material,

    fn hit(sphere: Sphere, r: Ray, tMin: f64, tMax: f64, record: *HitRecord) bool {
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
        if (root < tMin or tMax < root) {
            root = (-half_b + sqrtd) / a;
            if (root < tMin or tMax < root) {
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

    fn boudingBox(sphere: Sphere, time0: f64, time1: f64, outputBox: *Aabb) bool {
        _ = time0;
        _ = time1;
        const o = sphere.center;
        const r = sphere.radius;
        outputBox.* = .{
            .min = o.sub(.{ .x = r, .y = r, .z = r }),
            .max = o.add(.{ .x = r, .y = r, .z = r }),
        };
        return true;
    }

    fn deinit(sphere: *const Sphere, allocator: anytype) void {
        allocator.destroy(sphere.material);
    }
};

const MovingSphere = struct {
    center0: Point3,
    center1: Point3,
    time0: f64,
    time1: f64,
    radius: f64,
    material: *const Material,

    fn hit(sphere: MovingSphere, r: Ray, tMin: f64, tMax: f64, record: *HitRecord) bool {
        const center_ = sphere.center(r.time);
        const oc = r.origin.sub(center_);
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
        if (root < tMin or tMax < root) {
            root = (-half_b + sqrtd) / a;
            if (root < tMin or tMax < root) {
                // out of range
                return false;
            }
        }

        record.t = root;
        record.p = r.at(root);
        const outward_normal = (record.p.sub(center_)).div(sphere.radius);
        record.front_face = Vec3.dot(outward_normal, r.dir) < 0.0;
        if (record.front_face) {
            record.normal = outward_normal;
        } else {
            record.normal = outward_normal.mul(-1.0);
        }
        record.material = sphere.material;

        return true;
    }

    fn boudingBox(sphere: MovingSphere, time0: f64, time1: f64, outputBox: *Aabb) bool {
        const o0 = sphere.center(time0);
        const o1 = sphere.center(time1);
        const r = sphere.radius;
        const box0 = .{
            .min = o0.sub(.{ .x = r, .y = r, .z = r }),
            .max = o0.add(.{ .x = r, .y = r, .z = r }),
        };
        const box1 = .{
            .min = o1.sub(.{ .x = r, .y = r, .z = r }),
            .max = o1.add(.{ .x = r, .y = r, .z = r }),
        };
        outputBox.* = surroundingBox(box0, box1);
        return true;
    }

    fn center(sphere: MovingSphere, t: f64) Point3 {
        return sphere.center0.add(sphere.center1.sub(sphere.center0).mul((t - sphere.time0) / (sphere.time1 - sphere.time0)));
    }

    fn deinit(sphere: *const MovingSphere, allocator: anytype) void {
        allocator.destroy(sphere.material);
    }
};

fn makeSphere(center: Point3, radius: f64, material: *const Material) Hittable {
    return .{
        .sphere = .{
            .center = center,
            .radius = radius,
            .material = material,
        },
    };
}

const HittableList = struct {
    objects: ArrayList(Hittable),

    fn hit(list: HittableList, r: Ray, tMin: f64, tMax: f64, record: *HitRecord) bool {
        var hit_anything = false;
        var closest_so_far = tMax;

        for (list.objects.items) |object| {
            var rec: HitRecord = undefined;
            if (object.hit(r, tMin, closest_so_far, &rec)) {
                hit_anything = true;
                closest_so_far = rec.t;
                record.* = rec;
            }
        }
        return hit_anything;
    }

    fn boudingBox(list: HittableList, time0: f64, time1: f64, outputBox: *Aabb) bool {
        if (list.objects.items.len == 0) {
            return false;
        }
        var firstBox = true;
        var tmpBox: Aabb = undefined;
        for (list.objects.items) |object| {
            if (!object.boudingBox(time0, time1, tmpBox)) {
                return false;
            }
            outputBox = if (firstBox) tmpBox else surroundingBox(outputBox, tmpBox);
            firstBox = false;
        }
        return true;
    }

    fn deinit(list: *const HittableList, allocator: anytype) void {
        for (list.objects.items) |object| {
            object.deinit(allocator);
        }
        list.objects.deinit();
    }
};

const inf = math.inf(f64);
const pi = math.pi;

fn deg2rad(degree: f64) f64 {
    return degree * pi / 180.0;
}

// [min, max)
fn randomInt(comptime T: type, rand: std.rand.Random, min: T, max: T) T {
    return rand.intRangeLessThan(T, min, max);
}

// [0, 1)
fn randomReal01(rand: std.rand.Random) f64 {
    return rand.float(f64);
}

// [min, max)
fn randomReal(rand: std.rand.Random, min: f64, max: f64) f64 {
    return min + randomReal01(rand) * (max - min);
}

const Aabb = struct {
    min: Point3,
    max: Point3,

    fn hit(aabb: Aabb, r: Ray, tMin: f64, tMax: f64) bool {
        var tMin_ = tMin;
        var tMax_ = tMax;
        {
            const s0 = (aabb.min.x - r.origin.x) / r.dir.x;
            const s1 = (aabb.max.x - r.origin.x) / r.dir.x;
            const t0 = @min(s0, s1);
            const t1 = @max(s0, s1);
            tMin_ = @max(t0, tMin_);
            tMax_ = @min(t1, tMax_);
            if (tMax_ <= tMin_) {
                return false;
            }
        }
        {
            const s0 = (aabb.min.y - r.origin.y) / r.dir.y;
            const s1 = (aabb.max.y - r.origin.y) / r.dir.y;
            const t0 = @min(s0, s1);
            const t1 = @max(s0, s1);
            tMin_ = @max(t0, tMin_);
            tMax_ = @min(t1, tMax_);
            if (tMax_ <= tMin_) {
                return false;
            }
        }
        {
            const s0 = (aabb.min.z - r.origin.z) / r.dir.z;
            const s1 = (aabb.max.z - r.origin.z) / r.dir.z;
            const t0 = @min(s0, s1);
            const t1 = @max(s0, s1);
            tMin_ = @max(t0, tMin_);
            tMax_ = @min(t1, tMax_);
            if (tMax_ <= tMin_) {
                return false;
            }
        }
        return true;
    }
};

fn surroundingBox(box0: Aabb, box1: Aabb) Aabb {
    return .{
        .min = .{
            .x = @min(box0.min.x, box1.min.x),
            .y = @min(box0.min.y, box1.min.y),
            .z = @min(box0.min.z, box1.min.z),
        },
        .max = .{
            .x = @max(box0.max.x, box1.max.x),
            .y = @max(box0.max.y, box1.max.y),
            .z = @max(box0.max.z, box1.max.z),
        },
    };
}

const BvhNode = struct {
    left: *Hittable, // TODO
    right: *Hittable, // TODO
    box: Aabb,

    fn init(
        objects: ArrayList(*Hittable),
        start: usize,
        end: usize,
        time0: f64,
        time1: f64,
    ) BvhNode {
        var left: *Hittable = undefined;
        var right: *Hittable = undefined;

        var objects_ = objects;
        const axis = randomInt(u8, 0, 3);

        const objectSpan = end - start;
        if (objectSpan == 1) {
            left = objects.items[start];
            right = objects.items[start];
        } else if (objectSpan == 2) {
            if (BvhNode.boxCompare(axis, objects[start], objects[start + 1])) {
                left = objects[start];
                right = objects[start + 1];
            } else {
                left = objects[start + 1];
                right = objects[start];
            }
        } else {
            std.sort.sort(*Hittable, objects_, axis, BvhNode.boxCompare);
            const mid = start + objectSpan / 2;
            left.* = .{ .bvhNode = BvhNode.init(objects, start, mid, time0, time1) };
            right.* = .{ .bvhNode = BvhNode.init(objects, mid, end, time0, time1) };
        }

        var boxLeft: Aabb = undefined;
        var boxRight: Aabb = undefined;

        if (!left.boudingBox(time0, time1, boxLeft) || !right.boudingBox(time0, time1, boxRight)) {
            // ERROR
        }

        return .{
            .left = left,
            .right = right,
            .box = surroundingBox(boxLeft, boxRight),
        };
    }

    fn boxCompare(axis: u8, a: *const Hittable, b: *const Hittable) bool {
        if (axis == 0) {
            return a.x < b.x;
        } else if (axis == 1) {
            return a.y < b.y;
        } else {
            return a.z < b.z;
        }
    }

    fn hit(node: BvhNode, r: Ray, tMin: f64, tMax: f64, record: *HitRecord) bool {
        if (!node.box.hit(r, tMin, tMax)) {
            return false;
        }

        const hitLeft = node.left.hit(r, tMin, tMax, record);
        const hitRight = node.right.hit(r, tMin, if (hitLeft) record.t else tMax, record);

        return hitLeft or hitRight;
    }

    fn boudingBox(node: BvhNode, time0: f64, time1: f64, outputBox: *Aabb) bool {
        _ = time0;
        _ = time1;
        outputBox.* = node.box;
        return true;
    }

    fn deinit(node: *const BvhNode, allocator: anytype) void {
        _ = node;
        _ = allocator;
    }
};

const TextureTag = enum {
    solid,
};

const Texture = union(TextureTag) {
    solid: SolidTexture,

    fn value(tx: Texture, u: f64, v: f64, p: Vec3) Color {
        return switch (tx) {
            TextureTag.solid => |solidTx| solidTx.value(u, v, p),
        };
    }
};

const SolidTexture = struct {
    color: Color,

    fn value(tx: SolidTexture, u: f64, v: f64, p: Vec3) Color {
        _ = u;
        _ = v;
        _ = p;
        return tx.color;
    }
};

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

    fn getRay(camera: Camera, rand: std.rand.Random, s: f64, t: f64) Ray {
        const rd = randomPointInUnitDisk(rand).mul(camera.lens_radius);
        const offset = camera.u.mul(rd.x).add(camera.v.mul(rd.y));
        const dir = camera.lower_left_corner.add(camera.horizontal.mul(s)).add(camera.vertical.mul(t)).sub(camera.origin).sub(offset);
        return .{
            .origin = camera.origin.add(offset),
            .dir = dir,
            .time = randomReal(rand, camera.time0, camera.time1),
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

fn generateRandomScene(rand: std.rand.Random, allocator: anytype) !Hittable {
    var hittable_objects = ArrayList(Hittable).init(allocator);

    var mat_ground = try allocator.create(Material);
    var mat1 = try allocator.create(Material);
    var mat2 = try allocator.create(Material);
    var mat3 = try allocator.create(Material);

    mat_ground.* = .{ .diffuse = .{ .albedo = rgb(0.5, 0.5, 0.5) } };
    mat1.* = .{ .dielectric = .{ .ir = 1.5 } };
    mat2.* = .{ .diffuse = .{ .albedo = rgb(0.4, 0.2, 0.1) } };
    mat3.* = .{ .metal = .{ .albedo = rgb(0.7, 0.6, 0.5), .fuzz = 0.0 } };

    try hittable_objects.append(makeSphere(.{ .x = 0, .y = -1000, .z = 0 }, 1000, mat_ground));
    try hittable_objects.append(makeSphere(.{ .x = 0, .y = 1, .z = 0 }, 1.0, mat1));
    try hittable_objects.append(makeSphere(.{ .x = -4, .y = 1, .z = 0 }, 1.0, mat2));
    try hittable_objects.append(makeSphere(.{ .x = 4, .y = 1, .z = 0 }, 1.0, mat3));

    var a: i32 = -3;
    while (a < 3) : (a += 1) {
        var b: i32 = -3;
        while (b < 3) : (b += 1) {
            const choose_mat = randomReal01(rand);
            const center = Point3{
                .x = @intToFloat(f64, a) + 0.9 * randomReal01(rand),
                .y = 0.2,
                .z = @intToFloat(f64, b) + 0.9 * randomReal01(rand),
            };

            if (center.sub(.{ .x = 4, .y = 0.2, .z = 0 }).norm() <= 0.9) {
                continue;
            }

            var mat_sphere = try allocator.create(Material);
            if (choose_mat < 0.8) {
                // diffuse
                const albedo = Color.random01(rand).mulV(Color.random01(rand));
                mat_sphere.* = .{ .diffuse = .{ .albedo = albedo } };
                const center1 = center.add(.{ .x = 0, .y = randomReal(rand, 0, 0.5), .z = 0 });
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
                const albedo = Color.random(rand, 0.5, 1);
                const fuzz = randomReal(rand, 0, 0.5);
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

    var rng = std.rand.DefaultPrng.init(42);
    var rand = rng.random();

    // Image
    const aspect_ratio = 3.0 / 2.0;
    const image_width = 600;
    const image_height = @floatToInt(comptime_int, @divTrunc(image_width, aspect_ratio));
    const samples_per_pixel = 50;
    const max_depth = 50;

    // World
    const world = try generateRandomScene(rand, allocator);
    defer world.deinit(allocator);

    // Camera
    const camera = Camera.init(
        .{ .x = 12, .y = 2, .z = 3 },
        .{ .x = 0, .y = 0, .z = 0 },
        .{ .x = 0, .y = 1, .z = 0 },
        20.0,
        aspect_ratio,
        0.1,
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
