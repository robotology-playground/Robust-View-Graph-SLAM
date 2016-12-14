function C = pwg_to_C(pwg)
% With simple constraint output, do this first:
%   Cn = apply_nominal_R(C);
%   Cnn = constraint_loop_statistics(Cn, 1);

k = 1;
G = zeros(max(size(pwg)));
for i = 1:size(pwg,1)
    for j = 1:size(pwg,2)
        if isempty(pwg(i,j).z), continue, end
        G(i,j) = 1;
        C(k).edge = [i j];
        C(k).z = pwg(i,j).z;
        C(k).R = pwg(i,j).R;
        C(k).w = pwg(i,j).yes;
        k = k+1;
    end
end