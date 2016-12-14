function J = numerical_jacobian_i(model, normalise, idx, offset, varargin)
%function J = numerical_jacobian_i(model, normalise, i, offset, ...)
%
% INPUTS:
%   model - handle to function, y = model(x1, x2, ..., xi, ...)
%   normalise - handle to normalisation function, e = normalise(y1 - y2)
%   i - index of model argument xi about which to compute dy/dxi
%   offset - (optional) finite difference for approximating hyperplane
%   ... - arguments for model (x1, ..., xi, ...)
%
% OUTPUT:
%   J - Jacobian dy/dxi computed about the point (x1, x2, ...)
%
% REMARKS:
%   1. Function handle arguments may be passed as @model or 'model'. 
%   2. The normalise handle is optional; use only for discontinuous models. 
%   3. The offset argument is optional (defaults to 1e-9).
%   3. Jacobian matrix is computed via central differencing. This is a
%   quick and nasty way to compute derivatives and one should consider
%   automatic differentiation if more accurate results are required.
%   4. This function is essentially the same as numerical_jacobian_cd.m but
%   has an improved interface. The older versions are depreciated.
%
% Tim Bailey 2009.

if nargin == 3 || isempty(offset), offset = 1e-9; end

x = varargin{idx};
y = feval(model, varargin{:});
lenx = length(x);
leny = length(y);
J = zeros(leny, lenx);

for i=1:lenx
    xu = x(i) + offset;
    xl = x(i) - offset;

    varargin{idx}(i) = xu;
    yu = feval(model, varargin{:}); 
    varargin{idx}(i) = xl;
    yl = feval(model, varargin{:});    
    varargin{idx}(i) = x(i);
    
    dy = yu - yl;
    if ~isempty(normalise)
        dy = feval(normalise, dy);
    end
    
    J(:,i) = dy/(xu - xl);  
    % Numerically better to divide by xu - xl rather than 2*offset. 
end
