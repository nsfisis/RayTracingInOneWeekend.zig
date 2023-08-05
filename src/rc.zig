const std = @import("std");

pub fn Rc(comptime T: type) type {
    const Cell = struct {
        value: T,
        ref_count: usize,
    };

    return struct {
        const Self = @This();

        allocator: std.mem.Allocator,
        ptr: *Cell,

        pub fn init(allocator: std.mem.Allocator) !Self {
            const ptr = try allocator.create(Cell);
            ptr.ref_count = 1;
            return .{
                .allocator = allocator,
                .ptr = ptr,
            };
        }

        pub fn deinit(self: *const Self) void {
            self.ptr.ref_count -= 1;
            if (self.ptr.ref_count == 0) {
                self.allocator.destroy(self.ptr);
            }
        }

        pub fn clone(self: *const Self) Self {
            self.ptr.ref_count += 1;
            return .{
                .allocator = self.allocator,
                .ptr = self.ptr,
            };
        }

        pub fn get(self: *const Self) *const T {
            return &self.ptr.value;
        }

        pub fn get_mut(self: *Self) *T {
            return &self.ptr.value;
        }
    };
}
