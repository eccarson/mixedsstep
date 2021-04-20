% basisparams.m
% This file computes basis parameters for specified polynomials.
%
% Input:
%   s : the dimension of the basis to be computed
%   A : matrix for which we will compute basis coefficients
%   basis_info: a struct which has the quantity 'type' defined, which is a
%   string specifying the polynomials; options are 'monomial', 'newton', or
%   'chebyshev'. If this string is anything else, a monomial basis is
%   specified.
%
% Output:
%   alp,bet,gam: Vectors which store the computed coefficients
%   T: tridiagonal matrix containins alp, bet, gam on diagonals
%
% Last edited by: Erin Carson, 2021
%

function [alp,bet,gam, T] = basisparams(s, A, basis_info)


if(strcmp(basis_info.type, 'newton'))
    
    % Obtain extremal eigenvalues of A. Note that by commenting out 'eig' and
    % uncommenting the lines that use 'eigs', the basis can be generated using only
    % estimates of extremal eigenvalues (although need to set an appropriate
    % number of iterations). Has not been thoroughly tested.
    ee = eig(full(A));
    mx = max(ee);
    mn = min(ee);
    % mx = eigs(A,1,'LM');
    % mn = eigs(A,1,'SM');
    
    % Don't use scaled newton polynomials; edit to change this.
    basis_info.scale = ones(s,1);
    
    % Use max and min ritz values to compute leja points
    bbb = lejapoints(s, mn, mx);
    
    alp = bbb;
    bet = zeros(s-1,1);
    gam = ones(s,1);
    
elseif(strcmp(basis_info.type, 'chebyshev'))
    
    % Obtain extremal eigenvalues of A. Note that by commenting out 'eig' and
    % uncommenting the lines that use 'eigs', the basis can be generated using only
    % estimates of extremal eigenvalues (although need to set an appropriate
    % number of iterations). Has not been thoroughly tested.
    ee = eig(full(A));
    mx = max(ee);
    mn = min(ee);
    % mx = eigs(A,1,'LM');
    % mn = eigs(A,1,'SM');
    
    cc =(mx+mn)/2;
    aa = abs(mx-cc);
    bb = 0;
    dd = sqrt(aa^2-bb^2);
    gamm = max(aa,bb);
        
    alp = cc.*ones(s,1);
    bet = (dd)^2/(4*gamm).*ones(s-1,1);
    gam = [2*gamm; gamm.*ones(s-1,1)];
    
    
else % Assume monomial basis
    alp = zeros(s,1);
    bet = zeros(s-1,1);
    gam = ones(s,1);
end


T = diag(alp,0) + diag(bet,1) + diag(gam(1:end-1),-1);
T = [T;zeros(1,s-1),gam(end)];



