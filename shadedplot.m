%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Name:    shadeplot.m        Created: 11/06/16    Revised: 11/14/16
%
%% Usage:   Generate figure 4,5 in the multiscale network paper
%
%% Inputs:  coef_vs_mat := coefficient matrix
%           tt    := title text
%
%% Output:  
%% Calls:   Only internal Matlab functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = shadedplot(coef_vs_mat, tt)

shadedErrorBar(ceil((1:1499)*5/3)-500, mean(abs(coef_vs_mat(1,1,:,:)), 3),std(abs(coef_vs_mat(1,1,:,:)),0,3)/5)
title(tt);
xlabel('Time(ms)');
ylabel({'$\|\theta\|_2$'}, 'Interpreter','latex');
