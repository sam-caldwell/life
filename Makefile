# Simple wrapper around CMake

BUILD_DIR ?= build
CONFIG ?= Release

.PHONY: build clean configure

build: configure
	cmake --build $(BUILD_DIR) --config $(CONFIG) -j

configure:
	cmake -S . -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=$(CONFIG)

clean:
	rm -rf $(BUILD_DIR)
	mkdir -p $(BUILD_DIR)

