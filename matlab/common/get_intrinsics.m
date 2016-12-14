% decide which intrinsics
function [K1, kc1, K2, kc2] = get_intrinsics(options, i, j)
if nargin < 2; error('At least two inputs (options and camera i), are required'); end
if isodd(i)
    K1 = options.K1;
    if nargout > 1
        kc1 = options.kc1;
    end
else
    K1 = options.K2;
    if nargout > 1
        kc1 = options.kc2;
    end
end
if nargin < 3
    return
end
if isodd(j)
    K2 = options.K1;
    if nargout > 3
        kc2 = options.kc1;
    end
else
    K2 = options.K2;
    if nargout > 3
        kc2 = options.kc2;
    end
end

function b = isodd(i)
b = false;
if mod(i,2)
    b = true;
end

%function b = iseven(i)
%if mod(i,2)
%    b = false;
%else
%    b = even;
%end