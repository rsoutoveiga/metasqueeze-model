all: debug release

debug:
	@echo "config debug type"
	cmake -S . -B build/debug -D CMAKE_BUILD_TYPE=Debug
	@echo "build debug type"
	cmake --build build/debug --config Debug

release:
	@echo "config release type"
	cmake -S . -B build/release -D CMAKE_BUILD_TYPE=Release
	@echo "build release type"
	cmake --build build/release --config Release

build_debug:
	@echo "build debug type"
	cmake --build build/debug --config Debug

build_release:
	@echo "build release type"
	cmake --build build/release --config Release

clean:
	cmake --build build/debug/ --target clean
	cmake --build build/release/ --target clean

clean_debug:
	cmake --build build/debug/ --target clean

clean_release:
	cmake --build build/release/ --target clean


