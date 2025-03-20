# lajolla
UCSD CSE 272 renderer

The primary goal of our project is to enhance the LaJolla renderer by integrating Gradient-Domain
Path Tracing [1] algorithm into its regular path tracing implementation. This integration will build
upon the existing Monte Carlo sampling-based path tracing by incorporating gradient-domain sam-
pling, shift mapping and poisson reconstruction. By combining gradient information with tradition-
ally sampled pixels, gradient-domain sampling preserves the low-frequency details of conventional
images while leveraging the lower variance in sampled gradients to reduce high-frequency noise.
Additionally, we aim to replicate some of the authors’ analyses to compare renders produced by
gradient-domain path tracing and traditional path tracing to validate the correctness of our imple-
mentation.

Main files where we made changes :

1. src/render.cpp : Contains the poisson solver function \& our gradient\_path\_render function that performs GDPT.
2. src/path\_tracing.h : Contains our implementation of simultaneous evaluation of base and 4 offset paths.
3. src/intersection.h : Custom struct defined to hold each path contribution, probablities and multiple importance sampling weights.
4. src/parsers/parse\_scene.cpp : Updated the scene parser to map "gradpath" option to Integrator::GradPath object.
5. src/scene.h : Added GradPath to Integrator enum.
6. Also added some renders in cbox\_gdpt, cbox\_path, custom\_gdpt\_final\_render.exr and gdpt\_renders.




# Build
All the dependencies are included. Use CMake to build.
If you are on Unix systems, try
```
mkdir build
cd build
cmake ..
cmake --build .
```
It requires compilers that support C++17 (gcc version >= 8, clang version >= 7, Apple Clang version >= 11.0, MSVC version >= 19.14).

Apple M1 users: you might need to build Embree from scratch since the prebuilt MacOS binary provided is built for x86 machines. (But try build command above first.)

# Run
Try 
```
cd build
./lajolla ../scenes/cbox/cbox.xml
```
This will generate an image "image.exr".

To view the image, use [hdrview](https://github.com/wkjarosz/hdrview), or [tev](https://github.com/Tom94/tev).

# Acknowledgement
The renderer is heavily inspired by [pbrt](https://pbr-book.org/), [mitsuba](http://www.mitsuba-renderer.org/index_old.html), and [SmallVCM](http://www.smallvcm.com/).

We use [Embree](https://www.embree.org/) for ray casting.

We use [pugixml](https://pugixml.org/) to parse XML files.

We use [pcg](https://www.pcg-random.org/) for random number generation.

We use [stb_image](https://github.com/nothings/stb) and [tinyexr](https://github.com/syoyo/tinyexr) for reading & writing images.

We use [miniz](https://github.com/richgel999/miniz) for compression & decompression.

We use [tinyply](https://github.com/ddiakopoulos/tinyply) for parsing PLY files.

Many scenes in the scenes folder are directly downloaded from [http://www.mitsuba-renderer.org/download.html](http://www.mitsuba-renderer.org/download.html). Scenes courtesy of Wenzel Jakob, Cornell Program of Computer Graphics, Marko Dabrovic, Eric Veach, Jonas Pilo, and Bernhard Vogl.
