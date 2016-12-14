function heading_o = head_mod(heading)

% head_mod.m - heading angle modulation function

if heading > pi
    heading_o = heading - 2*pi;
elseif heading < -pi
    heading_o = heading + 2*pi;
else
    heading_o = heading;
end