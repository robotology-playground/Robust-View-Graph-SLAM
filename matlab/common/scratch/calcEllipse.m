function [X,Y] = calcEllipse(varargin)
% function [X,Y] = calculateEllipse(x, y, a, b, angle, steps)
%# This functions returns points to draw an ellipse
%#
%# @param x X coordinate
%# @param y Y coordinate
%# @param a Semimajor axis
%# @param b Semiminor axis
%# @param angle Angle of the ellipse (in rad)
%#
% Source: http://stackoverflow.com/questions/2153768/draw-ellipse-and-ellipsoid-in-matlab/24531259#24531259
% Modified by Christian Fässler

steps = 360;

if nargin == 1 || nargin == 2
x = varargin{1}.X0_in;
y = varargin{1}.Y0_in;
a = varargin{1}.a;
b = varargin{1}.b;
angle = varargin{1}.phi;
if nargin == 2
steps = varargin{2};
end
else if nargin == 5 || nargin == 6
x = varargin{1};
y = varargin{2};
a = varargin{3};
b = varargin{4};
angle = varargin{5};
if nargin == 6
steps = varargin{6};
end
else
error('Wrong input');
end
end

beta = -angle;
sinbeta = sin(beta);
cosbeta = cos(beta);

alpha = linspace(0, 2*pi, steps)';
sinalpha = sin(alpha);
cosalpha = cos(alpha);

X = round(x + (a * cosalpha * cosbeta - b * sinalpha * sinbeta));
Y = round(y + (a * cosalpha * sinbeta + b * sinalpha * cosbeta));

if nargout==1, X = [X Y]; end
end