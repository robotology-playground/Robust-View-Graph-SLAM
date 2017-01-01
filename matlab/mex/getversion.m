%--------------------------------------------------------------------------
function v = getversion
% determine the MATLAB version, and return it as a double.
v = sscanf(version, '%d.%d.%d');
v = 10.^(0:-1:-(length(v)-1)) * v;