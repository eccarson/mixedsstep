# mixedsstep
MATLAB implementations of mixed precision s-step Lanczos and CG, using double and simulated quadruple precision.

## Related publications
* Carson, E., & Gergelits, T. (2021). [Mixed Precision $s$-step Lanczos and Conjugate Gradient Algorithms](https://arxiv.org/pdf/2103.09210.pdf). arXiv preprint arXiv:2103.09210.

## Included MATLAB files
* **_cg.m_** is a simple implementation of the Hestenes and Stiefel variant of the classical conjugate gradient algorithm
* **_cacg.m_** is an implementation of the s-step (also called communication avoiding) variant of the CG method, using uniform precision (double)
* **_cacg_mixed.m_** is a mixed precision implementation of the s-step CG algorithm, which uses quadruple precision for computing and applying the Gram matrix and double precision elsewhere
* **_ca_lanczos_bounds.m_** is an implementation of the s-step (also called communication avoiding) variant of the Lanczos method, using uniform precision (double). This code also monitors the size of various quantities and also computes quantities needed to demonstrate the error bounds in the associated paper. 
* **_mixed_ca_lanczos_bounds.m_** is a mixed precision implementation of the s-step Lanczos algorithm, which uses quadruple precision for computing and applying the Gram matrix and double precision elsewhere. This code also monitors the size of various quantities and also computes quantities needed to demonstrate the error bounds in the associated paper. 
* **_computeBasis.m_** computes an s-dimensional Krylov subspace basis using the specified matrix, vector, and polynomial coefficients
* **_basisparams.m_** computes the polynomial basis coefficients used in the computeBasis function (options are monomial, Newton, or Chebyshev polynomials)
* **_lejapoints.m_** computes a set of leja points of specified size within a specified interval; used for computing the Newton basis coefficients
* **_compplots.m_** generates plots showing measured quantities and bounds for uniform precision and mixed precision s-step Lanczos for a given problem
* **_compplotscg.m_** generates convergence plots (in terms of relative error in the A-norm) for cg, cacg, and cacg_mixed for a given problem
* **_strakosmatrix.m_** generates a diagonal matrix with specified parameters (dimension, largest and smallest eigenvalue, and a parameter that controls clustering of the spectrum), commonly used in testing numerical properties of CG/Lanczos
* **_lanczos_example.m_** is an example driver for running Lanczos tests using this code. In particular, this script runs tests to compare uniform and mixed precision s-step Lanczos, using the diagonal test problem, s=6, and a monomial basis. This generates the plots that appear in Figure 5.1 of the associated preprint. 
* **_cg_example.m_** is an example driver for running CG tests using this code.  In particular, this script runs tests to compare the convergence of classical CG and uniform and mixed precision s-step CG, using the diagonal test problem, s=5, and a monomial basis. This generates the upper-right subplot in Figure 7.1 of the associated preprint. 
* The **_vpa/_** directory contains equivalent versions of these files that work with MATLAB's built-in vpa instead of the Advanpix toolkit. 

## Requirements
* The codes have been developed and tested with MATLAB 2020a.
* Double precision (built-in MATLAB type) is used as the working precision, and we use the Advanpix Multiprecision Computing Toolbox for extended precision computations. A free trial of Advanpix is available for download from https://www.advanpix.com/.
* If you prefer not to use the Advanpix library, use the code within the vpa directory, which uses MATLAB's vpa to perform higher precision computations. Note that this is much slower than Advanpix.  

For examples of how to use this code, please see the drivers cg_example.m and lanczos_example.m.

## License
See license.txt for licensing information.
