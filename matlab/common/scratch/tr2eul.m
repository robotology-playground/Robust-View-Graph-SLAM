%TR2EUL Convert a homogeneous transform matrix to Euler angle form
%
%	[PHI THETA PSI] = TR2EUL(M)
%
% Returns a vector of roll/pitch/yaw angles corresponding to M, either a rotation
% matrix or the rotation part of a homogeneous transform.
% The 3 angles correspond to rotations about the Z, Y and Z axes respectively.
%
% See also:  EUL2TR, TR2RPY

% Copyright (C) 1993-2008, by Peter I. Corke
%
% This file is part of The Robotics Toolbox for Matlab (RTB).
%
% RTB is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% RTB is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
%
% You should have received a copy of the GNU Leser General Public License
% along with RTB.  If not, see <http://www.gnu.org/licenses/>.

function euler = tr2eul(m)

s = size(m);
if length(s) > 2,
    euler = [];
    for i=1:s(3),
        euler = [euler; tr2eul(m(:,:,i))];
    end
    return
end

euler = zeros(1,3);

% Method as per Paul, p 69.
% phi = atan2(ay, ax)
% Only positive phi is returned.
if abs(m(1,3)) < eps & abs(m(2,3)) < eps,
    % singularity
    euler(1) = 0;
    sp = 0;
    cp = 1;
    euler(2) = atan2(cp*m(1,3) + sp*m(2,3), m(3,3));
    euler(3) = atan2(-sp * m(1,1) + cp * m(2,1), -sp*m(1,2) + cp*m(2,2));
else
    euler(1) = atan2(m(2,3), m(1,3));
    sp = sin(euler(1));
    cp = cos(euler(1));
    euler(2) = atan2(cp*m(1,3) + sp*m(2,3), m(3,3));
    euler(3) = atan2(-sp * m(1,1) + cp * m(2,1), -sp*m(1,2) + cp*m(2,2));
end
