% strakosmatrix.m
% Generate diagonal "strakos matrix" test matrix
%
% Input:
%   n : dimension of matrix
%   l1 and ln: are smallest and largest eigenvalues, resp.
%   p: 0<p<=1 determines eigenvalue distribution
%
% Last edited by: Erin Carson, 2021
%

function A = strakosmatrix(n, l1, ln, p)

d(1) = l1;
d(n) = ln;
for i = 2:n-1
    d(i) = l1 + ((i-1)/(n-1))*(ln-l1)*p^(n-i);
end
A = diag(d);