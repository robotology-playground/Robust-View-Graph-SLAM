function CG = C_to_CG(C)
% converts constraints C into camera graph CG
N = max(max(vertcat(C.edge)));
CG = zeros(N);
for i = 1:length(C);
    CG(C(i).edge(1), C(i).edge(2)) = 1;
end