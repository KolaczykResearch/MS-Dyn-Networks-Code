%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Name:    plotbp.m        Created: 11/06/16    Revised: 11/14/16
%
%% Usage:   Generate figure 3 in the multiscale network paper
%
%% Inputs:  store_bp := detected change points
%           color    := color used for the bar
%
%% Output:  
%% Calls:   Only internal Matlab functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h] = plotbp(store_bp, color)

if nargin < 2
    color = [1 .5 0];
end
 
bp = ceil(cell2mat(store_bp)*5/3) -500;

h = histogram(bp, 150); % in version after 2014b, use function histogram instead of hist
h.FaceColor = color; 

line([0 0], [0 20], 'Color', [0 0 0], 'LineStyle','--');
line([1000 1000], [0 20], 'Color', [0 0 0],'LineStyle','--');