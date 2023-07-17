const Ray = @import("ray.zig").Ray;
const Point3 = @import("vec.zig").Point3;

pub const Aabb = struct {
    min: Point3,
    max: Point3,

    pub fn hit(aabb: Aabb, r: Ray, t_min: f64, t_max: f64) bool {
        var t_min_ = t_min;
        var t_max_ = t_max;
        {
            const s0 = (aabb.min.x - r.origin.x) / r.dir.x;
            const s1 = (aabb.max.x - r.origin.x) / r.dir.x;
            const t0 = @min(s0, s1);
            const t1 = @max(s0, s1);
            t_min_ = @max(t0, t_min_);
            t_max_ = @min(t1, t_max_);
            if (t_max_ <= t_min_) {
                return false;
            }
        }
        {
            const s0 = (aabb.min.y - r.origin.y) / r.dir.y;
            const s1 = (aabb.max.y - r.origin.y) / r.dir.y;
            const t0 = @min(s0, s1);
            const t1 = @max(s0, s1);
            t_min_ = @max(t0, t_min_);
            t_max_ = @min(t1, t_max_);
            if (t_max_ <= t_min_) {
                return false;
            }
        }
        {
            const s0 = (aabb.min.z - r.origin.z) / r.dir.z;
            const s1 = (aabb.max.z - r.origin.z) / r.dir.z;
            const t0 = @min(s0, s1);
            const t1 = @max(s0, s1);
            t_min_ = @max(t0, t_min_);
            t_max_ = @min(t1, t_max_);
            if (t_max_ <= t_min_) {
                return false;
            }
        }
        return true;
    }

    pub fn surroundingBox(box0: Aabb, box1: Aabb) Aabb {
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
};
