# Image Denoiser

This project contains MATLAB scripts implementing various numerical methods for image denoising. The algorithms aim to reduce noise in images using iterative solvers and matrix operations.

## Files Description

- **CG.m**  
  Implements the Conjugate Gradient method for solving linear systems related to image denoising.

- **Denoise.m**  
  Main script to apply the denoising algorithms to input images.

- **FormMatrix.m**  
  Constructs the system matrix used in the numerical methods.

- **FormRHS.m**  
  Forms the right-hand side vector for the linear system.

- **GS.m**  
  Implements the Gauss-Seidel iterative solver.

- **Jacobi.m**  
  Implements the Jacobi iterative solver.

- **SOR.m**  
  Implements the Successive Over-Relaxation method.

- **set_image.m**  
  Script for setting up the input image and parameters.

- **test.m**  
  Testing script to validate and compare the performance of different solvers.

- **README.md**  
  This file.

## How to Use

1. Set the input image and parameters in `set_image.m`.
2. Run `Denoise.m` to apply the denoising process.
3. Use `test.m` to test and compare different methods.

## Requirements

- MATLAB or GNU Octave

*Created by Rukhan4*
