# Gradient-Based-Path-Tracing
Implementation of the paper Gradient-domain path tracing in LaJolla (CSE272)

The primary goal of our project is to enhance the LaJolla renderer by integrating Gradient-Domain
Path Tracing algorithm into its regular path tracing implementation. This integration will build
upon the existing Monte Carlo sampling-based path tracing by incorporating gradient-domain sam-
pling, shift mapping and poisson reconstruction. By combining gradient information with tradition-
ally sampled pixels, gradient-domain sampling preserves the low-frequency details of conventional
images while leveraging the lower variance in sampled gradients to reduce high-frequency noise.
Additionally, we aim to replicate some of the authors’ analyses to compare renders produced by
gradient-domain path tracing and traditional path tracing to validate the correctness of our imple-
mentation
