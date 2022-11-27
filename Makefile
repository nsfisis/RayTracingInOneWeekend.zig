.PHONY: all
all: fmt
	@rm -f out.ppm
	@zig build run > out.ppm

.PHONY: fmt
fmt:
	@zig fmt .
