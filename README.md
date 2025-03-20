# Gradient-Based-Path-Tracing
Implementation of the paper Gradient-domain path tracing in LaJolla (CSE272)

The primary goal of our project is to enhance the LaJolla renderer by integrating Gradient-Domain
Path Tracing algorithm into its regular path tracing implementation. This integration will build
upon the existing Monte Carlo sampling-based path tracing by incorporating gradient-domain sampling, 
shift mapping and poisson reconstruction. By combining gradient information with traditionally sampled pixels,
gradient-domain sampling preserves the low-frequency details of conventional
images while leveraging the lower variance in sampled gradients to reduce high-frequency noise.
Additionally, we aim to replicate some of the authors’ analyses to compare renders produced by
gradient-domain path tracing and traditional path tracing to validate the correctness of our implementation


Main files where we made changes :

1. src/render.cpp : Contains the poisson solver function & our gradient_path_render function that performs GDPT.
2. src/path_tracing.h : Contains our implementation of simultaneous evaluation of base and 4 offset paths.
3. src/intersection.h : Custom struct defined to hold each path contribution, probablities and multiple importance sampling weights.
4. src/parsers/parse_scene.cpp : Updated the scene parser to map "gradpath" option to Integrator::GradPath object.
5. src/scene.h : Added GradPath to Integrator enum.
6. Also added some renders in cbox_gdpt, cbox_path, custom_gdpt_final_render.exr and gdpt_renders.


