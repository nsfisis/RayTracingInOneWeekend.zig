const std = @import("std");
const debug = std.debug;
const math = std.math;
const pi = math.pi;
const ArrayList = std.ArrayList;

const Rc = @import("../rc.zig").Rc;

const Ray = @import("ray.zig").Ray;
const Vec3 = @import("vec.zig").Vec3;
const Point3 = @import("vec.zig").Point3;
const Color = @import("vec.zig").Color;
const rgb = @import("vec.zig").rgb;
const randomPointInUnitSphere = @import("rand.zig").randomPointInUnitSphere;
const randomPointInUnitDisk = @import("rand.zig").randomPointInUnitDisk;
const randomUnitVector = @import("rand.zig").randomUnitVector;
const randomInt = @import("rand.zig").randomInt;
const randomReal01 = @import("rand.zig").randomReal01;
const randomReal = @import("rand.zig").randomReal;

const Material = @import("material.zig").Material;

const HitRecord = @import("hit_record.zig").HitRecord;
const Aabb = @import("aabb.zig").Aabb;

const HittableTag = enum {
    sphere,
    movingSphere,
    list,
    bvhNode,
    xyRect,
    xzRect,
    yzRect,
    box,
    translate,
    rotateY,
};

pub const Hittable = union(HittableTag) {
    const Self = @This();

    sphere: Sphere,
    movingSphere: MovingSphere,
    list: HittableList,
    bvhNode: BvhNode,
    xyRect: XyRect,
    xzRect: XzRect,
    yzRect: YzRect,
    box: Box,
    translate: Translate,
    rotateY: RotateY,

    pub fn makeBox(p0: Point3, p1: Point3, material: Rc(Material), allocator: anytype) !Self {
        return .{ .box = try Box.init(p0, p1, material, allocator) };
    }

    pub fn translate(obj: Rc(Self), offset: Vec3) Self {
        return .{ .translate = .{ .object = obj, .offset = offset } };
    }

    pub fn rotateY(obj: Rc(Self), angle: f64) Self {
        return .{ .rotateY = RotateY.init(obj, angle) };
    }

    pub fn hit(self: Self, r: Ray, t_min: f64, t_max: f64, record: *HitRecord) bool {
        return switch (self) {
            HittableTag.sphere => |slf| slf.hit(r, t_min, t_max, record),
            HittableTag.movingSphere => |slf| slf.hit(r, t_min, t_max, record),
            HittableTag.list => |slf| slf.hit(r, t_min, t_max, record),
            HittableTag.bvhNode => |slf| slf.hit(r, t_min, t_max, record),
            HittableTag.xyRect => |slf| slf.hit(r, t_min, t_max, record),
            HittableTag.xzRect => |slf| slf.hit(r, t_min, t_max, record),
            HittableTag.yzRect => |slf| slf.hit(r, t_min, t_max, record),
            HittableTag.box => |slf| slf.hit(r, t_min, t_max, record),
            HittableTag.translate => |slf| slf.hit(r, t_min, t_max, record),
            HittableTag.rotateY => |slf| slf.hit(r, t_min, t_max, record),
        };
    }

    fn boudingBox(self: Self, time0: f64, time1: f64, output_box: *Aabb) bool {
        return switch (self) {
            HittableTag.sphere => |slf| slf.boudingBox(time0, time1, output_box),
            HittableTag.movingSphere => |slf| slf.boudingBox(time0, time1, output_box),
            HittableTag.list => |slf| slf.boudingBox(time0, time1, output_box),
            HittableTag.bvhNode => |slf| slf.boudingBox(time0, time1, output_box),
            HittableTag.xyRect => |slf| slf.boudingBox(time0, time1, output_box),
            HittableTag.xzRect => |slf| slf.boudingBox(time0, time1, output_box),
            HittableTag.yzRect => |slf| slf.boudingBox(time0, time1, output_box),
            HittableTag.box => |slf| slf.boudingBox(time0, time1, output_box),
            HittableTag.translate => |slf| slf.boudingBox(time0, time1, output_box),
            HittableTag.rotateY => |slf| slf.boudingBox(time0, time1, output_box),
        };
    }

    pub fn deinit(self: *const Self) void {
        return switch (self.*) {
            HittableTag.sphere => |slf| slf.deinit(),
            HittableTag.movingSphere => |slf| slf.deinit(),
            HittableTag.list => |slf| slf.deinit(),
            HittableTag.bvhNode => |slf| slf.deinit(),
            HittableTag.xyRect => |slf| slf.deinit(),
            HittableTag.xzRect => |slf| slf.deinit(),
            HittableTag.yzRect => |slf| slf.deinit(),
            HittableTag.box => |slf| slf.deinit(),
            HittableTag.translate => |slf| slf.deinit(),
            HittableTag.rotateY => |slf| slf.deinit(),
        };
    }
};

const Sphere = struct {
    center: Point3,
    radius: f64,
    material: Rc(Material),

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
        Sphere.getSphereUv(outward_normal, &record.u, &record.v);
        record.material = sphere.material.get();

        return true;
    }

    fn boudingBox(sphere: Sphere, time0: f64, time1: f64, output_box: *Aabb) bool {
        _ = time0;
        _ = time1;
        const o = sphere.center;
        const r = sphere.radius;
        output_box.* = .{
            .min = o.sub(.{ .x = r, .y = r, .z = r }),
            .max = o.add(.{ .x = r, .y = r, .z = r }),
        };
        return true;
    }

    fn getSphereUv(p: Point3, u: *f64, v: *f64) void {
        const phi = math.atan2(-p.z, p.x) + pi;
        const theta = math.acos(-p.y);
        u.* = phi / (2.0 * pi);
        v.* = theta / pi;
    }

    fn deinit(sphere: *const Sphere) void {
        sphere.material.deinit();
    }
};

const MovingSphere = struct {
    center0: Point3,
    center1: Point3,
    time0: f64,
    time1: f64,
    radius: f64,
    material: Rc(Material),

    fn hit(sphere: MovingSphere, r: Ray, t_min: f64, t_max: f64, record: *HitRecord) bool {
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
        if (root < t_min or t_max < root) {
            root = (-half_b + sqrtd) / a;
            if (root < t_min or t_max < root) {
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
        record.material = sphere.material.get();

        return true;
    }

    fn boudingBox(sphere: MovingSphere, time0: f64, time1: f64, output_box: *Aabb) bool {
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
        output_box.* = Aabb.surroundingBox(box0, box1);
        return true;
    }

    fn center(sphere: MovingSphere, t: f64) Point3 {
        return sphere.center0.add(sphere.center1.sub(sphere.center0).mul((t - sphere.time0) / (sphere.time1 - sphere.time0)));
    }

    fn deinit(sphere: *const MovingSphere) void {
        sphere.material.deinit();
    }
};

const HittableList = struct {
    objects: ArrayList(Hittable),

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

    fn boudingBox(list: HittableList, time0: f64, time1: f64, output_box: *Aabb) bool {
        if (list.objects.items.len == 0) {
            return false;
        }
        var first_box = true;
        var tmp_box: Aabb = undefined;
        for (list.objects.items) |object| {
            if (!object.boudingBox(time0, time1, &tmp_box)) {
                return false;
            }
            output_box.* = if (first_box) tmp_box else Aabb.surroundingBox(output_box.*, tmp_box);
            first_box = false;
        }
        return true;
    }

    fn deinit(list: *const HittableList) void {
        for (list.objects.items) |object| {
            object.deinit();
        }
        list.objects.deinit();
    }
};

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

        const objects_ = objects;
        const axis = randomInt(u8, 0, 3);

        const object_span = end - start;
        if (object_span == 1) {
            left = objects.items[start];
            right = objects.items[start];
        } else if (object_span == 2) {
            if (BvhNode.boxCompare(axis, objects[start], objects[start + 1])) {
                left = objects[start];
                right = objects[start + 1];
            } else {
                left = objects[start + 1];
                right = objects[start];
            }
        } else {
            std.sort.sort(*Hittable, objects_, axis, BvhNode.boxCompare);
            const mid = start + object_span / 2;
            left.* = .{ .bvhNode = BvhNode.init(objects, start, mid, time0, time1) };
            right.* = .{ .bvhNode = BvhNode.init(objects, mid, end, time0, time1) };
        }

        const box_left: Aabb = undefined;
        const box_right: Aabb = undefined;

        if (!left.boudingBox(time0, time1, box_left) || !right.boudingBox(time0, time1, box_right)) {
            // ERROR
        }

        return .{
            .left = left,
            .right = right,
            .box = Aabb.surroundingBox(box_left, box_right),
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

    fn hit(node: BvhNode, r: Ray, t_min: f64, t_max: f64, record: *HitRecord) bool {
        if (!node.box.hit(r, t_min, t_max)) {
            return false;
        }

        const hit_left = node.left.hit(r, t_min, t_max, record);
        const hit_right = node.right.hit(r, t_min, if (hit_left) record.t else t_max, record);

        return hit_left or hit_right;
    }

    fn boudingBox(node: BvhNode, time0: f64, time1: f64, output_box: *Aabb) bool {
        _ = time0;
        _ = time1;
        output_box.* = node.box;
        return true;
    }

    fn deinit(node: *const BvhNode) void {
        node.left.deinit();
        node.right.deinit();
    }
};

const XyRect = struct {
    x0: f64,
    x1: f64,
    y0: f64,
    y1: f64,
    k: f64,
    material: Rc(Material),

    pub fn hit(rect: XyRect, r: Ray, t_min: f64, t_max: f64, record: *HitRecord) bool {
        const t = (rect.k - r.origin.z) / r.dir.z;
        if (t < t_min or t > t_max) {
            return false;
        }
        const x = r.origin.x + t * r.dir.x;
        const y = r.origin.y + t * r.dir.y;
        if (x < rect.x0 or x > rect.x1 or y < rect.y0 or y > rect.y1) {
            return false;
        }

        record.u = (x - rect.x0) / (rect.x1 - rect.x0);
        record.v = (y - rect.y0) / (rect.y1 - rect.y0);
        record.t = t;
        record.material = rect.material.get();
        record.p = r.at(t);

        const outward_normal = Vec3{ .x = 0, .y = 0, .z = 1 };
        record.front_face = Vec3.dot(outward_normal, r.dir) < 0.0;
        if (record.front_face) {
            record.normal = outward_normal;
        } else {
            record.normal = outward_normal.mul(-1.0);
        }
        return true;
    }

    pub fn boudingBox(rect: XyRect, time0: f64, time1: f64, output_box: *Aabb) bool {
        _ = time0;
        _ = time1;

        // The bounding box must have non-zero width in each dimension, so pad the Z
        // dimension a small amount.
        output_box.* = .{
            .min = .{ .x = rect.x0, .y = rect.y0, .z = rect.k - 0.0001 },
            .max = .{ .x = rect.x1, .y = rect.y1, .z = rect.k + 0.0001 },
        };
        return true;
    }

    pub fn deinit(rect: *const XyRect) void {
        rect.material.deinit();
    }
};

const XzRect = struct {
    x0: f64,
    x1: f64,
    z0: f64,
    z1: f64,
    k: f64,
    material: Rc(Material),

    pub fn hit(rect: XzRect, r: Ray, t_min: f64, t_max: f64, record: *HitRecord) bool {
        const t = (rect.k - r.origin.y) / r.dir.y;
        if (t < t_min or t > t_max) {
            return false;
        }
        const x = r.origin.x + t * r.dir.x;
        const z = r.origin.z + t * r.dir.z;
        if (x < rect.x0 or x > rect.x1 or z < rect.z0 or z > rect.z1) {
            return false;
        }

        record.u = (x - rect.x0) / (rect.x1 - rect.x0);
        record.v = (z - rect.z0) / (rect.z1 - rect.z0);
        record.t = t;
        record.material = rect.material.get();
        record.p = r.at(t);

        const outward_normal = Vec3{ .x = 0, .y = 1, .z = 0 };
        record.front_face = Vec3.dot(outward_normal, r.dir) < 0.0;
        if (record.front_face) {
            record.normal = outward_normal;
        } else {
            record.normal = outward_normal.mul(-1.0);
        }
        return true;
    }

    pub fn boudingBox(rect: XzRect, time0: f64, time1: f64, output_box: *Aabb) bool {
        _ = time0;
        _ = time1;

        // The bounding box must have non-zero width in each dimension, so pad the Y
        // dimension a small amount.
        output_box.* = .{
            .min = .{ .x = rect.x0, .y = rect.k - 0.0001, .z = rect.z0 },
            .max = .{ .x = rect.x1, .y = rect.k + 0.0001, .z = rect.z1 },
        };
        return true;
    }

    pub fn deinit(rect: *const XzRect) void {
        rect.material.deinit();
    }
};

const YzRect = struct {
    y0: f64,
    y1: f64,
    z0: f64,
    z1: f64,
    k: f64,
    material: Rc(Material),

    pub fn hit(rect: YzRect, r: Ray, t_min: f64, t_max: f64, record: *HitRecord) bool {
        const t = (rect.k - r.origin.x) / r.dir.x;
        if (t < t_min or t > t_max) {
            return false;
        }
        const y = r.origin.y + t * r.dir.y;
        const z = r.origin.z + t * r.dir.z;
        if (y < rect.y0 or y > rect.y1 or z < rect.z0 or z > rect.z1) {
            return false;
        }

        record.u = (y - rect.y0) / (rect.y1 - rect.y0);
        record.v = (z - rect.z0) / (rect.z1 - rect.z0);
        record.t = t;
        record.material = rect.material.get();
        record.p = r.at(t);

        const outward_normal = Vec3{ .x = 1, .y = 0, .z = 0 };
        record.front_face = Vec3.dot(outward_normal, r.dir) < 0.0;
        if (record.front_face) {
            record.normal = outward_normal;
        } else {
            record.normal = outward_normal.mul(-1.0);
        }
        return true;
    }

    pub fn boudingBox(rect: YzRect, time0: f64, time1: f64, output_box: *Aabb) bool {
        _ = time0;
        _ = time1;

        // The bounding box must have non-zero width in each dimension, so pad the X
        // dimension a small amount.
        output_box.* = .{
            .min = .{ .x = rect.k - 0.0001, .y = rect.y0, .z = rect.z0 },
            .max = .{ .x = rect.k + 0.0001, .y = rect.y1, .z = rect.z1 },
        };
        return true;
    }

    pub fn deinit(rect: *const YzRect) void {
        rect.material.deinit();
    }
};

const Box = struct {
    min: Point3,
    max: Point3,
    sides: HittableList,

    pub fn init(p0: Point3, p1: Point3, material: Rc(Material), allocator: anytype) !Box {
        var sides = ArrayList(Hittable).init(allocator);

        try sides.append(.{ .xyRect = .{ .x0 = p0.x, .x1 = p1.x, .y0 = p0.y, .y1 = p1.y, .k = p1.z, .material = material } });
        try sides.append(.{ .xyRect = .{ .x0 = p0.x, .x1 = p1.x, .y0 = p0.y, .y1 = p1.y, .k = p0.z, .material = material.clone() } });
        try sides.append(.{ .xzRect = .{ .x0 = p0.x, .x1 = p1.x, .z0 = p0.z, .z1 = p1.z, .k = p1.y, .material = material.clone() } });
        try sides.append(.{ .xzRect = .{ .x0 = p0.x, .x1 = p1.x, .z0 = p0.z, .z1 = p1.z, .k = p0.y, .material = material.clone() } });
        try sides.append(.{ .yzRect = .{ .y0 = p0.y, .y1 = p1.y, .z0 = p0.z, .z1 = p1.z, .k = p1.x, .material = material.clone() } });
        try sides.append(.{ .yzRect = .{ .y0 = p0.y, .y1 = p1.y, .z0 = p0.z, .z1 = p1.z, .k = p0.x, .material = material.clone() } });

        return .{
            .min = p0,
            .max = p1,
            .sides = .{
                .objects = sides,
            },
        };
    }

    pub fn hit(box: Box, r: Ray, t_min: f64, t_max: f64, record: *HitRecord) bool {
        return box.sides.hit(r, t_min, t_max, record);
    }

    pub fn boudingBox(box: Box, time0: f64, time1: f64, output_box: *Aabb) bool {
        _ = time0;
        _ = time1;
        output_box.* = .{
            .min = box.min,
            .max = box.max,
        };
        return true;
    }

    pub fn deinit(box: *const Box) void {
        box.sides.deinit();
    }
};

const Translate = struct {
    const Self = @This();

    object: Rc(Hittable),
    offset: Vec3,

    pub fn hit(self: Self, r: Ray, t_min: f64, t_max: f64, record: *HitRecord) bool {
        const offset_r: Ray = .{
            .origin = r.origin.sub(self.offset),
            .dir = r.dir,
            .time = r.time,
        };
        if (!self.object.get().hit(offset_r, t_min, t_max, record)) {
            return false;
        }
        record.p = record.p.add(self.offset);
        return true;
    }

    pub fn boudingBox(self: Self, time0: f64, time1: f64, output_box: *Aabb) bool {
        const result = self.object.get().boudingBox(time0, time1, output_box);
        output_box.* = .{
            .min = output_box.min.add(self.offset),
            .max = output_box.max.add(self.offset),
        };
        return result;
    }

    pub fn deinit(self: *const Self) void {
        self.object.deinit();
    }
};

const RotateY = struct {
    const Self = @This();

    object: Rc(Hittable),
    sin_t: f64,
    cos_t: f64,
    bbox: Aabb,

    pub fn init(object: Rc(Hittable), t: f64) Self {
        const sin_t = math.sin(t);
        const cos_t = math.cos(t);
        var bbox: Aabb = undefined;
        _ = object.get().boudingBox(0, 0, &bbox);

        var min: Point3 = .{ .x = math.inf(f64), .y = math.inf(f64), .z = math.inf(f64) };
        var max: Point3 = .{ .x = -math.inf(f64), .y = -math.inf(f64), .z = -math.inf(f64) };

        var i: u32 = 0;
        while (i < 2) : (i += 1) {
            var j: u32 = 0;
            while (j < 2) : (j += 1) {
                var k: u32 = 0;
                while (k < 2) : (k += 1) {
                    const i_f: f64 = @floatFromInt(i);
                    const j_f: f64 = @floatFromInt(j);
                    const k_f: f64 = @floatFromInt(k);
                    const x = i_f * bbox.max.x + (1.0 - i_f) * bbox.min.x;
                    const y = j_f * bbox.max.y + (1.0 - j_f) * bbox.min.y;
                    const z = k_f * bbox.max.z + (1.0 - k_f) * bbox.min.z;

                    const nx = cos_t * x + sin_t * z;
                    const nz = -sin_t * x + cos_t * z;

                    const tester: Vec3 = .{ .x = nx, .y = y, .z = nz };

                    min.x = @min(min.x, tester.x);
                    min.y = @min(min.y, tester.y);
                    min.z = @min(min.z, tester.z);
                    max.x = @max(max.x, tester.x);
                    max.y = @max(max.y, tester.y);
                    max.z = @max(max.z, tester.z);
                }
            }
        }

        return .{
            .object = object,
            .sin_t = sin_t,
            .cos_t = cos_t,
            .bbox = .{ .min = min, .max = max },
        };
    }

    pub fn hit(self: Self, r: Ray, t_min: f64, t_max: f64, record: *HitRecord) bool {
        // Change the ray from world space to object space.
        var origin = r.origin;
        var dir = r.dir;

        origin.x = self.cos_t * r.origin.x - self.sin_t * r.origin.z;
        origin.z = self.sin_t * r.origin.x + self.cos_t * r.origin.z;

        dir.x = self.cos_t * r.dir.x - self.sin_t * r.dir.z;
        dir.z = self.sin_t * r.dir.x + self.cos_t * r.dir.z;

        const rotated_r: Ray = .{
            .origin = origin,
            .dir = dir,
            .time = r.time,
        };

        if (!self.object.get().hit(rotated_r, t_min, t_max, record)) {
            return false;
        }

        // Change the intersection point from object space to world space.
        // Rotate -t around axis Y.
        //   cos(-t) = cos(t)
        //   sin(-t) = -sin(t)
        var p = record.p;
        p.x = self.cos_t * record.p.x + self.sin_t * record.p.z;
        p.z = -self.sin_t * record.p.x + self.cos_t * record.p.z;

        // Change the normal from object space to world space.
        var normal = record.normal;
        normal.x = self.cos_t * record.normal.x + self.sin_t * record.normal.z;
        normal.z = -self.sin_t * record.normal.x + self.cos_t * record.normal.z;

        record.p = p;
        record.normal = normal;

        return true;
    }

    pub fn boudingBox(self: Self, time0: f64, time1: f64, output_box: *Aabb) bool {
        _ = time0;
        _ = time1;
        output_box.* = self.bbox;
        return true;
    }

    pub fn deinit(self: *const Self) void {
        self.object.deinit();
    }
};
