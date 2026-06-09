# SPHinXsys installation on Apple Silicon macOS

This instruction sheet records a working installation path for SPHinXsys on an M1/M2/M3 Apple Silicon Mac, using CMake 3.30.1, Apple Clang 17, Ninja, ccache, Homebrew gfortran, and vcpkg `2024.11.16`.

It follows the official SPHinXsys macOS flow, with one additional workaround for a `libffi` build failure seen with newer Apple Clang / macOS SDK combinations.

## Tested environment

The successful build path was verified on:

```bash
git version 2.50.1 Apple Git-155
Python 3.13.9
cmake version 3.30.1
Apple clang version 17.0.0
Target: arm64-apple-darwin25.3.0
GNU Fortran Homebrew GCC 15.2.0_1
ninja 1.13.2
ccache 4.13.4
```

The machine was native Apple Silicon:

```bash
uname -m
# arm64
```

## Important rules

1. Use `arm64-osx`, not `x64-osx`.
2. Use CMake 3.x. SPHinXsys requires CMake `>= 3.16` and `< 4.0`.
3. Use vcpkg branch/tag `2024.11.16`, as requested by SPHinXsys.
4. Do not reuse stale `build/` directories after a failed CMake configure.
5. Do not run the SYCL / oneAPI / CUDA section on macOS. That section is for Linux GPU builds.

## Recommended workspace layout

This guide uses:

```bash
/Users/adheshsagar/CodeFiles/bfm11/
├── vcpkg
└── sphinxsys
```

For another user, replace:

```bash
/Users/adheshsagar/CodeFiles/bfm11
```

with their own workspace path.

Set the workspace variable:

```bash
export WORK="$HOME/CodeFiles/bfm11"
mkdir -p "$WORK"
cd "$WORK"
```

## 1. Install system dependencies

Install Xcode Command Line Tools first if needed:

```bash
xcode-select --install
```

Install Homebrew dependencies:

```bash
brew update

brew install \
  git \
  python \
  pkg-config \
  ccache \
  gcc \
  ninja \
  autoconf \
  automake \
  autoconf-archive
```

Check tools:

```bash
uname -m
git --version
python3 --version
cmake --version
clang --version
gfortran --version
ninja --version
ccache --version
```

Expected:

```bash
uname -m
# arm64

cmake --version
# cmake version 3.30.1, or another 3.x version before 4.0
```

## 2. Clone pinned vcpkg

```bash
cd "$WORK"

git clone -b 2024.11.16 https://github.com/microsoft/vcpkg.git
cd "$WORK/vcpkg"

./bootstrap-vcpkg.sh -disableMetrics
./vcpkg version
```

Expected vcpkg version observed:

```text
vcpkg package management program version 2024-11-12-eb492805e92a2c14a230f5c3deb3e89f6771c321
```

Set Apple Silicon triplets:

```bash
export VCPKG_ROOT="$WORK/vcpkg"
export VCPKG_DEFAULT_TRIPLET=arm64-osx
export VCPKG_DEFAULT_HOST_TRIPLET=arm64-osx
```

## 3. Install dependencies, first normal attempt

Run:

```bash
cd "$WORK/vcpkg"

./vcpkg install --clean-after-build --allow-unsupported \
  eigen3:arm64-osx \
  tbb:arm64-osx \
  boost-program-options:arm64-osx \
  boost-geometry:arm64-osx \
  simbody:arm64-osx \
  gtest:arm64-osx \
  pybind11:arm64-osx \
  spdlog:arm64-osx
```

If all packages install successfully, skip to section 5.

## 4. Apple Clang 17 `libffi` workaround

Use this section only if `libffi:arm64-osx` fails with an error like:

```text
src/aarch64/sysv.S
error: invalid CFI advance_loc expression
.cfi_def_cfa x1, 40;
.cfi_adjust_cfa_offset ...
error: building libffi:arm64-osx failed with: BUILD_FAILED
```

The workaround that succeeded was:

- Keep vcpkg pinned at `2024.11.16`.
- Override only the `libffi` port with newer vcpkg’s `libffi 3.5.2`.
- Add the newer helper ports `vcpkg-make` and `vcpkg-cmake-get-vars`.
- Patch the newer `libffi` portfile so the old pinned vcpkg installs into the expected package folders.

### 4.1 Create a repair branch

```bash
cd "$WORK/vcpkg"

git switch sphinxsys-libffi-newer-port 2>/dev/null || git switch -c sphinxsys-libffi-newer-port
git fetch origin master
```

### 4.2 Bring in newer port directories

```bash
git checkout origin/master -- ports/libffi
git checkout origin/master -- ports/vcpkg-make
git checkout origin/master -- ports/vcpkg-cmake-get-vars

git add ports/libffi ports/vcpkg-make ports/vcpkg-cmake-get-vars
git commit -m "Use newer libffi and helper ports for Apple Clang 17"
```

If Git says `nothing to commit`, continue.

Verify the newer libffi port:

```bash
grep -n '"version"' ports/libffi/vcpkg.json
```

Expected:

```text
"version": "3.5.2"
```

### 4.3 Patch `ports/libffi/portfile.cmake`

Run this whole block from the vcpkg root:

```bash
cd "$WORK/vcpkg"

python3 <<'PY'
from pathlib import Path

path = Path("ports/libffi/portfile.cmake")
s = path.read_text()

start = s.index("vcpkg_make_configure(")

depth = 0
end = None
for i in range(start, len(s)):
    if s[i] == "(":
        depth += 1
    elif s[i] == ")":
        depth -= 1
        if depth == 0:
            end = i
            break

if end is None:
    raise RuntimeError("Could not find vcpkg_make_configure block")

old_block = s[start:end+1]

# Remove any previous OPTIONS_RELEASE / OPTIONS_DEBUG section if present.
cut = old_block.find("    OPTIONS_RELEASE")
if cut != -1:
    base_block = old_block[:cut].rstrip()
else:
    base_block = old_block[:-1].rstrip()

new_block = base_block + r'''
    OPTIONS_RELEASE
        "--prefix=/"
        "--bindir=/tools/${PORT}/bin"
        "--sbindir=/tools/${PORT}/sbin"
        "--libdir=/lib"
        "--includedir=/include"
        "--datarootdir=/share/${PORT}"
    OPTIONS_DEBUG
        "--prefix=/debug"
        "--bindir=/tools/${PORT}/debug/bin"
        "--sbindir=/tools/${PORT}/debug/sbin"
        "--libdir=/debug/lib"
        "--includedir=/include"
        "--datarootdir=/debug/share/${PORT}"
)'''

s = s[:start] + new_block + s[end+1:]
path.write_text(s)

print("Repatched ports/libffi/portfile.cmake with DESTDIR-safe prefixes")
PY
```

Verify that the patch does not contain absolute package paths:

```bash
grep -n -- '--prefix' ports/libffi/portfile.cmake
grep -n -- '--libdir' ports/libffi/portfile.cmake
grep -n -- '--includedir' ports/libffi/portfile.cmake
```

Expected lines include:

```text
"--prefix=/"
"--prefix=/debug"
"--libdir=/lib"
"--libdir=/debug/lib"
"--includedir=/include"
```

They should not contain:

```text
/Users/.../vcpkg/packages
```

Commit the patch:

```bash
git add ports/libffi/portfile.cmake
git commit -m "Fix libffi DESTDIR-safe install prefixes"
```

### 4.4 Clean old failed libffi package

```bash
cd "$WORK/vcpkg"

rm -rf buildtrees/libffi
rm -rf packages/libffi_arm64-osx
rm -rf installed/arm64-osx/share/libffi
rm -f installed/arm64-osx/include/ffi.h
rm -f installed/arm64-osx/include/ffitarget.h

find installed/arm64-osx/lib -name 'libffi*' -delete 2>/dev/null || true
find installed/arm64-osx/debug/lib -name 'libffi*' -delete 2>/dev/null || true
```

### 4.5 Install libffi alone

```bash
export VCPKG_ROOT="$WORK/vcpkg"
export VCPKG_DEFAULT_TRIPLET=arm64-osx
export VCPKG_DEFAULT_HOST_TRIPLET=arm64-osx

./vcpkg install --clean-after-build --allow-unsupported libffi:arm64-osx
```

Verify:

```bash
./vcpkg list | grep libffi
```

Expected:

```text
libffi:arm64-osx 3.5.2
```

### 4.6 Continue installing SPHinXsys dependencies

```bash
./vcpkg install --clean-after-build --allow-unsupported \
  eigen3:arm64-osx \
  tbb:arm64-osx \
  boost-program-options:arm64-osx \
  boost-geometry:arm64-osx \
  simbody:arm64-osx \
  gtest:arm64-osx \
  pybind11:arm64-osx \
  spdlog:arm64-osx
```

If a download fails with:

```text
curl: (6) Could not resolve host: github.com
```

it is a DNS/network issue, not a build issue. Test and retry:

```bash
curl -I https://github.com
./vcpkg install --clean-after-build --allow-unsupported spdlog:arm64-osx
```

Verify all needed packages:

```bash
./vcpkg list | grep -E 'eigen3|tbb|boost-program-options|boost-geometry|simbody|gtest|pybind11|spdlog|libffi'
```

## 5. Clone SPHinXsys

```bash
cd "$WORK"

git clone https://github.com/Xiangyu-Hu/SPHinXsys.git sphinxsys
cd "$WORK/sphinxsys"
```

If already cloned:

```bash
cd "$WORK/sphinxsys"
git pull
```

## 6. Configure SPHinXsys

Remove any old failed build folder:

```bash
cd "$WORK/sphinxsys"
rm -rf build
```

Configure with Ninja, clang, vcpkg, and explicit ARM settings:

```bash
cmake -G Ninja \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_C_COMPILER=/usr/bin/clang \
  -D CMAKE_CXX_COMPILER=/usr/bin/clang++ \
  -D CMAKE_OSX_ARCHITECTURES=arm64 \
  -D CMAKE_TOOLCHAIN_FILE="$WORK/vcpkg/scripts/buildsystems/vcpkg.cmake" \
  -D VCPKG_TARGET_TRIPLET=arm64-osx \
  -D VCPKG_HOST_TRIPLET=arm64-osx \
  -D TBB_DIR="$WORK/vcpkg/installed/arm64-osx/share/tbb" \
  -D CMAKE_C_COMPILER_LAUNCHER=ccache \
  -D CMAKE_CXX_COMPILER_LAUNCHER=ccache \
  -D TEST_STATE_RECORDING=OFF \
  -S . \
  -B build
```

The explicit `TBB_DIR` is included because the successful run initially failed at:

```text
Could not find a package configuration file provided by "TBB"
```

If that happens, check:

```bash
cd "$WORK/vcpkg"
./vcpkg list | grep -i tbb
find installed/arm64-osx -iname '*TBBConfig.cmake' -o -iname '*tbb-config.cmake'
```

Expected:

```text
installed/arm64-osx/share/tbb/TBBConfig.cmake
```

## 7. Build SPHinXsys

```bash
cd "$WORK/sphinxsys"

cmake --build build --config Release --verbose
```

A warning like this is acceptable if the build finishes:

```text
ld: warning: ignoring duplicate libraries: '-ldl', '-lm'
```

## 8. Run tests

The full test suite should be run sequentially:

```bash
cd "$WORK/sphinxsys/build"

ctest -j 1 --output-on-failure
```

Sequential testing matters because SPHinXsys uses all available cores internally.

## 9. Run one example manually

Because this guide uses Ninja, do not use `make -j 7`.

For example, to build and run `test_2d_dambreak`:

```bash
cd "$WORK/sphinxsys"

cmake --build build --target test_2d_dambreak

cd build/tests/2d_examples/test_2d_dambreak/bin
./test_2d_dambreak
```

## 10. What not to do on Apple Silicon macOS

Do not use:

```bash
./vcpkg env --triplet=x64-osx
```

Do not use:

```bash
x64-osx
```

Do not continue with the SYCL / oneAPI / CUDA instructions on macOS. That section is for Linux GPU/SYCL builds.

Do not keep retrying CMake configuration after failures without deleting:

```bash
rm -rf "$WORK/sphinxsys/build"
```

## 11. Troubleshooting summary

### Problem: `libffi` fails with `invalid CFI advance_loc expression`

Use section 4.

### Problem: `vcpkg-make does not exist`

You copied the newer `libffi` port but not the helper port. Run:

```bash
git checkout origin/master -- ports/vcpkg-make
```

### Problem: `vcpkg_cmake_get_vars was passed extra arguments: ADDITIONAL_LANGUAGES ASM`

You need the newer helper port:

```bash
git checkout origin/master -- ports/vcpkg-cmake-get-vars
```

Then remove old helper packages and retry:

```bash
./vcpkg remove vcpkg-make:arm64-osx vcpkg-cmake-get-vars:arm64-osx --recurse

rm -rf buildtrees/vcpkg-cmake-get-vars
rm -rf buildtrees/vcpkg-make
rm -rf buildtrees/libffi

rm -rf packages/vcpkg-cmake-get-vars_arm64-osx
rm -rf packages/vcpkg-make_arm64-osx
rm -rf packages/libffi_arm64-osx

rm -rf installed/arm64-osx/share/vcpkg-cmake-get-vars
rm -rf installed/arm64-osx/share/vcpkg-make
rm -rf installed/arm64-osx/share/libffi
```

### Problem: doubled package path in libffi

Example:

```text
packages/libffi_arm64-osx/Users/.../packages/libffi_arm64-osx/...
```

Your `libffi` portfile used absolute package paths. Reapply section 4.3 using the `/`, `/debug`, `/lib`, and `/include` prefixes.

### Problem: `Could not find TBBConfig.cmake`

Make sure TBB is installed and pass `TBB_DIR` explicitly:

```bash
-D TBB_DIR="$WORK/vcpkg/installed/arm64-osx/share/tbb"
```

### Problem: GitHub download fails

If the error is:

```text
curl: (6) Could not resolve host: github.com
```

check network/DNS and retry; do not clean all vcpkg packages.

```bash
curl -I https://github.com
```

## Source notes

- Official SPHinXsys installation page: https://www.sphinxsys.org/html/installation.html
- SPHinXsys requires CMake 3.16 or later but before 4.0, Python3, Git, and vcpkg branch `2024.11.16`.
- Official macOS instructions say M1/later Macs should use `arm64-osx` instead of `x64-osx`.
- Current vcpkg has `libffi 3.5.2`, with helper dependencies `vcpkg-cmake-get-vars` and `vcpkg-make`.
