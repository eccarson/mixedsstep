% lanczos_example_vpa.m
% This file shows an example script that generates plots for loss of
% orthogonality, normality, column error, and deviation of AV and T for
% uniform precision s-step Lanczos and mixed precision s-step Lanczos.
% The code will plot both the measured quantity and the bound. 
%
% Last edited by: Erin Carson, 2021
%

% Set matrix A; here A is a diagonal test problem with clustered
% eigenvalues of dimension 100
A = strakosmatrix(100,1e-3,100,.65);

% Set starting vector v
v = ones(size(A,1),1);
v = v./norm(v);

% Set value of s to use
s = 5;

% Set the type of polynomial basis you would like to use; options are
% 'monomial', 'chebyshev', or 'newton'
basis_info.type = 'monomial';

% Set the number of iterations to run
options.xlim = 100;

% Set the name of the matrix; this is for naming the output files
options.name = 'strakos';


% Create 'figs' folder for output if it doesn't exist yet
if ~exist('figs', 'dir')
    mkdir('figs')
end

% Call compplots to run the algorithms and generate the plots 
compplots_vpa(A, v, s, basis_info, options)