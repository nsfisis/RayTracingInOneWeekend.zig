.PHONY: all
all: fmt
	@rm -f out.ppm
	@zig build run > out.ppm

.PHONY: build
build: fmt
	@zig build

.PHONY: fmt
fmt:
	@zig fmt .
