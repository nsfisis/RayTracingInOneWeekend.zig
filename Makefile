BUILDFLAGS := -Doptimize=ReleaseFast

.PHONY: all
all: fmt
	@rm -f out.png
	@zig build ${BUILDFLAGS} run

.PHONY: build
build: fmt
	@zig build ${BUILDFLAGS}

.PHONY: fmt
fmt:
	@zig fmt .
