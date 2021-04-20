% compplots_vpa.m
% This file runs classical CG in double precision, uniform precision s-step CG in double, 
% and mixed precision s-step CG double/quad and creates plots convergence in terms of the relative
% A-norm of the error
%
% The zero vector is used as the initial approximate solution for all
% algorithms.
%
% Input:
%   A : a square, symmetric positive definite coefficient matrix
%   b : the right-hand side in Ax=b
%   s : the parameter s in s-step (s>0)
%   tol : the convergence tolerance
%   basis_info : a struct containing a quantity 'type' which defines the
%   polynomial basis to be used. Options are 'monomial', 'newton', or
%   'chebyshev'
%   options : a struct containing a quantity 'xlim', which gives the number
%   of iterations to perform, a quantity 'name', which is a string used for
%   naming the output file, and the quantity 'truesol', which is a vector
%   containing the true solution, used for measuring the error
%
%   See cg_example_vpa.m for example usage.
%
% Last edited by: Erin Carson, 2021
%
function compplotscg_vpa(A, b, s, basis_info, options)

% Run the classical CG algorithm with the specifiec parameters
results = cg_vpa(A, b, zeros(size(A,1),1), options);

% Run the uniform precision s-step CG algorithm with the specified
% parameters
resultsca = cacg_vpa(A, b, s, zeros(size(A,1),1), basis_info, options);

% Run the mixed precision s-step CG algorithm with the specified parameters
resultscam = cacg_mixed_vpa(A, b, s, zeros(size(A,1),1), basis_info, options);


% Set x axis limit
m = options.xlim;

% Generate and save convergence plot
h = figure()
semilogy(1:m-1, results.error_A_norm(1:m-1),'k-', 1:m-1, resultsca.error_A_norm(1:m-1),'r-', ...
    1:m-1, resultscam.error_A_norm(1:m-1),'b-');
tt = strcat('Relative Error, $s=$ ', num2str(s), ", ", basis_info.type, ' basis');
title(tt,'Interpreter','latex','FontSize',18);
legend('CG double', 'uniform $s$-step CG, double', 'mixed $s$-step CG, double/quad', 'Interpreter', 'latex','FontSize',14); 
xlabel('Iteration','Interpreter','latex','FontSize',14)
sname = strcat('figs/',options.name,'_cg_s',num2str(s),'_',basis_info.type);
saveas(h,sname,'fig')
saveas(h,sname,'pdf')
