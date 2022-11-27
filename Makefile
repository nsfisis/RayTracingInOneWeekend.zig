.PHONY: all
all:
	@zig build run

.PHONY: fmt
fmt:
	@zig fmt .
