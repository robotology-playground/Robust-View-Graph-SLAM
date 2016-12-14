function CG = convert_C_to_CG(C)
% With simple constraint output, do this first:
%   Cn = apply_nominal_R(C);
%   Cnn = constraint_loop_statistics(Cn, 1);

k = 1;
for i = 1:size(C,1)
    for j = 1:size(C,2)
        if isempty(C(i,j).z), continue, end
        CG(k).edge = [i j];
        CG(k).z = C(i,j).z;
        CG(k).R = C(i,j).R;
        CG(k).w = C(i,j).yes;
        k = k+1;
    end
end