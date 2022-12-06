BUILDFLAGS := -Drelease-fast=true

.PHONY: all
all: fmt
	@rm -f out.ppm
	@zig build ${BUILDFLAGS} run > out.ppm

.PHONY: build
build: fmt
	@zig build ${BUILDFLAGS}

.PHONY: fmt
fmt:
	@zig fmt .
