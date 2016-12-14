function p = line_plot_conversion_(s, t)
%function p = line_plot_conversion_(s, t)
%
% INPUTS:
%   s - (DxN) line start points
%   t - (DxN) line end points
%
% OUTPUT: 
%   p - points for D-dimensional line plotting
%
% Convert a list of lines so that they will be plotted as a set of
% unconnected lines but only require a single handle to do so. 
%
% Tim Bailey 2011.

[D,N] = size(s);
len = N*3;
p = zeros(D, len);
p(:, 1:3:end) = s;
p(:, 2:3:end) = t;
p(:, 3:3:end) = NaN;
