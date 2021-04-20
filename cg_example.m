% cg_example.m
% This file shows an example script that compares the convergence
% trajectories of classical CG in double precision, uniform precision
% s-step CG in double precision, and mixed precision s-step CG in
% double/quad precision. 
% The code will plot the relative error in the A-norm. 
%
% Last edited by: Erin Carson, 2021
%

% Set matrix A; here A is a diagonal test problem with clustered
% eigenvalues of dimension 100
A = strakosmatrix(100,1e-3,100,.65);

% Set the right-hand side b to be a unit vector with equal components in the eigenbasis of A
[V,D] = eig(full(A));
b = V*ones(size(A,1),1);
b = b./norm(b);

% Set value of s to use
s = 6;


% Set the type of polynomial basis you would like to use; options are
% 'monomial', 'chebyshev', or 'newton'
basis_info.type='monomial';

% Set the maximum number of iterations to run
options.xlim = 900;

% Set the name of the matrix; this is for naming the output files
options.name='strakos';

% Set the true solution by computing using backslash in quadruple precision
options.truesol = double(mp(A,34)\mp(b,34));

% Create folder for output if it doesn't exist yet
if ~exist('figs', 'dir')
    mkdir('figs')
end

% Call compplotscg to run the algorithms and generate the plots 
compplotscg(A, b, s, basis_info, options)
