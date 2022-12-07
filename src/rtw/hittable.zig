const std = @import("std");
const debug = std.debug;
const math = std.math;
const pi = math.pi;
const ArrayList = std.ArrayList;

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
const DiffuseMaterial = @import("material.zig").DiffuseMaterial;
const MetalMaterial = @import("material.zig").MetalMaterial;
const DielectricMaterial = @import("material.zig").DielectricMaterial;

const Texture = @import("texture.zig").Texture;
const SolidTexture = @import("texture.zig").SolidTexture;
const CheckerTexture = @import("texture.zig").CheckerTexture;

const HitRecord = @import("hit_record.zig").HitRecord;
const Aabb = @import("aabb.zig").Aabb;

const HittableTag = enum {
    sphere,
    movingSphere,
    list,
    bvhNode,
};

pub const Hittable = union(HittableTag) {
    sphere: Sphere,
    movingSphere: MovingSphere,
    list: HittableList,
    bvhNode: BvhNode,

    pub fn hit(h: Hittable, r: Ray, tMin: f64, tMax: f64, record: *HitRecord) bool {
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

    pub fn deinit(h: *const Hittable, allocator: anytype) void {
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
        Sphere.getSphereUv(outward_normal, &record.u, &record.v);
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

    fn getSphereUv(p: Point3, u: *f64, v: *f64) void {
        const phi = math.atan2(f64, -p.z, p.x) + pi;
        const theta = math.acos(-p.y);
        u.* = phi / (2.0 * pi);
        v.* = theta / pi;
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
        outputBox.* = Aabb.surroundingBox(box0, box1);
        return true;
    }

    fn center(sphere: MovingSphere, t: f64) Point3 {
        return sphere.center0.add(sphere.center1.sub(sphere.center0).mul((t - sphere.time0) / (sphere.time1 - sphere.time0)));
    }

    fn deinit(sphere: *const MovingSphere, allocator: anytype) void {
        allocator.destroy(sphere.material);
    }
};

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
            outputBox = if (firstBox) tmpBox else Aabb.surroundingBox(outputBox, tmpBox);
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
            .box = Aabb.surroundingBox(boxLeft, boxRight),
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
