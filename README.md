# mixedsstep
MATLAB implementations of mixed precision s-step Lanczos and CG

This repository contains MATLAB implementations of mixed precision s-step Lanczos and mixed precision s-step CG algorithms. This code was used to generate the plots in 
Carson, E., & Gergelits, T. (2021). Mixed Precision $s$-step Lanczos and Conjugate Gradient Algorithms. arXiv preprint arXiv:2103.09210.

Here, double precision is used as the working precision, and for extended precision, we use quadruple precision simulated via the Advanpix Multiprecision Toolbox. Running this code therefore requires that you have the Advanpix Toolbox installed (https://www.advanpix.com/). 

For examples of how to use this code, please see the drivers cg_example.m and lanczos_example.m.
