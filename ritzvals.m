% ritzvals.m
% This file shows an example script that generates plots showing the
% distance of computed Ritz values to the nearest eigenvalues of A after
% 100 iterations, for uniform precision s-step Lanczos, mixed precision 
% s-step Lanczos, and classical Lanczos for comparison.
%
% Last edited by: Erin Carson, 2021



% Set matrix A; here A is a diagonal test problem with clustered
% eigenvalues of dimension 100
A = strakosmatrix(100,1e-3,100,.65);

% Store eigenvalues of A (just the diagonal of the matrix here)
ee = diag(A);

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

% Run uniform precision s-step Lanczos with the specified parameters
[results, basis_info] = ca_lanczos_bounds(A, v, s, basis_info, options);

% Run mixed precision s-step Lanczos with the specified parameters
[resultsm, basis_infom] = mixed_ca_lanczos_bounds(A, v, s,basis_info, options);

% Run classical Lanczos with the specified parameters
[resultsc] = lanczos_alg(A, v, options);

% Compute the Ritz values for each run using extra precision
ritz = mp(eig(results.T),64);
ritzm = mp(eig(resultsm.T),64);
ritzc = mp(eig(resultsc.T),64);

% Find the distance to the closest eigenvalue of A
for i = 1:numel(ritz)
    [closest, dist] = findclosest(ritz(i), ee);
    ritzdist(i) = dist;
end

for i = 1:numel(ritzm)
    [closest, dist] = findclosest(ritzm(i), ee);
    ritzmdist(i) = dist;
end

for i = 1:numel(ritzc)
    [closest, dist] = findclosest(ritzc(i), ee);
    ritzcdist(i) = dist;
end

% Plot results
h = figure()
set(0,'DefaultAxesTitleFontWeight','normal');
axes1 = axes;
hold(axes1,'on');

semilogy(1:100, ritzdist,'bo', 'MarkerSize',10)
hold on
semilogy(1:100, ritzmdist,'rx','MarkerSize',10)
hold on
semilogy(1:100, ritzcdist,'k+','MarkerSize',10)

tt = strcat('Distance to eigenvalues of A, s = '," ", num2str(s), ", ", basis_info.type, ' basis');
title(tt,'Interpreter','tex','FontSize',16);
box(axes1,'on');
hold(axes1,'off');
axis([0 100, 1e-18, 10])
set(axes1,'FontSize',14,'YMinorTick','on','YScale','log');
sname = strcat('figs/',options.name,'_ritzvals',num2str(s),'_',basis_info.type);
saveas(h,sname,'fig')
saveas(h,sname,'pdf')
