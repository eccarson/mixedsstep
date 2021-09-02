% lanczos_slg.m
% This file runs the classical Lanczos algorithm (version which uses 
% two-term recurrences).
%
% Input:
%   A : a square, symmetric coefficient matrix
%   v : the starting vector for the Lanczos method
%   options : a struct containing a quantity 'xlim', which gives the number
%   of iterations to perform, and a quantity 'name', which is used for
%   naming the output file
%
%
% Last edited by: Erin Carson, 2021

function [results] = lanczos_alg(A, v, options)

% Get size of matrix and maximum number of nonzeros per row
n = size(A,1);
N = max(sum(A~=0,2));

% Set initial values for vectors (ensure unit starting vector)
v = v./norm(v);
v0 = v;
u0 = A*v0;
v(:,1) = v0;
u(:,1) = u0;

% Initialize quantities
beta(1) = 0;
m = 0;

% Begin the iterations!
while m < options.xlim    
        
        % Increment global iteration count
        m = m + 1;
                
        % Update iteration vectors
        alpha(m) = v(:,m)'*u(:,m);
        w = u(:,m) - alpha(m)*v(:,m);
        beta(m+1) = sqrt(w'*w);
        v(:,m+1) = w/beta(m+1);
        u(:,m+1) = A*v(:,m+1) -beta(m+1)*v(:,m);

end

% Store tridiagonal T matrix
results.T = diag(alpha(1:end),0)+diag(beta(2:end-1),1)+ diag(beta(2:end-1),-1);

end
