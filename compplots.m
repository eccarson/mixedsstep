% compplots.m
% This file runs uniform precision s-step Lanczos and mixed precision
% s-step Lanczos and creates plots of loss of
% orthogonality, normality, column error, and deviation of AV and T for the
% specified problem.
%
% Input:
%   A : a square, symmetric coefficient matrix
%   v : the starting vector for the Lanczos method
%   s : the parameter s in s-step (s>0)
%   basis_info : a struct containing a quantity 'type' which defines the
%   polynomial basis to be used. Options are 'monomial', 'newton', or
%   'chebyshev'
%   options : a struct containing a quantity 'xlim', which gives the number
%   of iterations to perform, and a quantity 'name', which is used for
%   naming the output file
%
%   See lanczos_example.m for example usage.
%
% Last edited by: Erin Carson, 2021
%
function compplots(A, v, s, basis_info, options)

% Run uniform precision s-step Lanczos with the specified parameters
[results, basis_info] = ca_lanczos_bounds(A, v, s, basis_info, options);

% Run mixed precision s-step Lanczos with the specified parameters
[resultsm, basis_infom] = mixed_ca_lanczos_bounds(A, v, s,basis_info, options);

% Set x axis limit
m = options.xlim+1;

%%%%%%%%%%%% Loss of orthogonality plot %%%%%%%%%%%
h = figure()
set(0,'DefaultAxesTitleFontWeight','normal');
 axes1 = axes;
 hold(axes1,'on');
semilogy(1:m-1, results.orthactual(1:m-1),'b-', 1:m-1, results.orthbound(1:m-1),'b--', ...
    1:m-1, resultsm.orthactual(1:m-1),'r-', 1:m-1, resultsm.orthbound(1:m-1),'r--', ...,
    1:m-1,results.Gamma_orth(1:m-1),'k-','LineWidth',1);
tt = strcat('Loss of orthogonality, s = ', " ",num2str(s), ", ", basis_info.type, ' basis');
title(tt,'Interpreter','none','FontSize',16);
% Uncomment below to plot legend
% legend('measured, uniform', 'bound, uniform', 'measured, mixed', 'bound, mixed','$\bar{\Gamma}_k$', 'Interpreter', 'latex','FontSize',14); 
xlabel('Iteration','Interpreter','none','FontSize',16)
box(axes1,'on');
hold(axes1,'off');
set(axes1,'FontSize',14,'YMinorTick','on','YScale','log');
sname = strcat('figs/',options.name,'_orthbound_s',num2str(s),'_',basis_info.type);
saveas(h,sname,'fig')
saveas(h,sname,'pdf')

%%%%%%%%%%%% Normality plot %%%%%%%%%%%
h = figure()
 axes1 = axes;
 hold(axes1,'on');
semilogy(1:m-1, results.normactual(1:m-1),'b-', 1:m-1, results.normbound(1:m-1),'b--', ...
    1:m-1, resultsm.normactual(1:m-1),'r-', 1:m-1, resultsm.normbound(1:m-1),'r--', ...
    1:m-1,results.Gamma_norm(1:m-1),'k-','LineWidth',1);
tt = strcat('Normality, s = ', " ",num2str(s), ", ", basis_info.type, ' basis');
title(tt,'Interpreter','none','FontSize',16);
set(0,'DefaultAxesTitleFontWeight','normal');
% Uncomment below to plot legend
 legend('measured, uniform', 'bound, uniform', 'measured, mixed', 'bound, mixed','$${\Gamma}_{sk+j}$$', 'Interpreter', 'none','FontSize',14); 
 legend1 = legend(axes1,'show');
set(legend1,'Interpreter','latex','FontSize',14);
 xlabel('Iteration','Interpreter','none','FontSize',16)
box(axes1,'on');
hold(axes1,'off');
 set(axes1,'FontSize',14,'YMinorTick','on','YScale','log');
sname = strcat('figs/',options.name,'_normbound_s',num2str(s),'_',basis_info.type);
saveas(h,sname,'fig')
saveas(h,sname,'pdf')

%%%%%%%%%%%% Column error plot %%%%%%%%%%%
h = figure()
 axes1 = axes;
 hold(axes1,'on');
semilogy(1:m-1, results.colactual(1:m-1),'b-', 1:m-1, results.colbound(1:m-1),'b--', ...
    1:m-1, resultsm.colactual(1:m-1),'r-', 1:m-1, resultsm.colbound(1:m-1),'r--', ...
    1:m-1,results.Gamma_colbound(1:m-1),'k-', 'LineWidth',1);
tt = strcat('Error in Columns, s = ', " ",num2str(s), ", ", basis_info.type, ' basis');
title(tt,'Interpreter','none','FontSize',16);
set(0,'DefaultAxesTitleFontWeight','normal');
% Uncomment below to plot legend
% legend('col, double', 'col bound, double', 'col, double/quad', 'col s-Lanc. double/quad','$\bar{\Gamma}_k$', 'Interpreter', 'latex','FontSize',14); 
xlabel('Iteration','Interpreter','none','FontSize',16)
box(axes1,'on');
hold(axes1,'off');
set(axes1,'FontSize',14,'YMinorTick','on','YScale','log');
sname = strcat('figs/',options.name,'_colbound_s',num2str(s),'_',basis_info.type);
saveas(h,sname,'fig')
saveas(h,sname,'pdf')

%%%%%%%%%%%% Deviation of AV and T plot %%%%%%%%%%%
h = figure()
 axes1 = axes;
 hold(axes1,'on');
semilogy(1:m-1, results.diffactual(1:m-1),'b-', 1:m-1, results.diffbound(1:m-1),'b--',...
    1:m-1, resultsm.diffactual(1:m-1),'r-', 1:m-1, resultsm.diffbound(1:m-1),'r--', ...
    1:m-1,results.Gamma_diffbound(1:m-1),'k-','LineWidth',1);
tt = strcat('Deviation of AV_m and T_m, s = '," ", num2str(s), ", ", basis_info.type, ' basis');
title(tt,'Interpreter','tex','FontSize',16);
set(0,'DefaultAxesTitleFontWeight','normal');
% Uncomment below to plot legend
% legend('diff, double', 'diff bound, double', 'diff, double/quad', 'diff s-Lanc. double/quad','$\bar{\Gamma}_k$', 'Interpreter', 'latex','FontSize',14); 
xlabel('Iteration','Interpreter','none','FontSize',16)
box(axes1,'on');
hold(axes1,'off');
set(axes1,'FontSize',14,'YMinorTick','on','YScale','log');
sname = strcat('figs/',options.name,'_diffbound_s',num2str(s),'_',basis_info.type);
saveas(h,sname,'fig')
saveas(h,sname,'pdf')
