function [pwg, G] = C_to_pwg(C)
edges = vertcat(C.edge);
N = max(edges(:));
%pwg = zeros(N);
G = zeros(N);
for i = 1:length(C);
    if G(C(i).edge(1),C(i).edge(2)) == 1, continue, end,
    if isfield(C, 'z')
        pwg(C(i).edge(1),C(i).edge(2)).z = C(i).z;
    else
        pwg(C(i).edge(1),C(i).edge(2)).z = [C(i).t;C(i).a];
    end
    if isfield(C, 'c');
        pwg(C(i).edge(1),C(i).edge(2)).c = C(i).c;
    elseif isfield(C, 'sw');
        pwg(C(i).edge(1),C(i).edge(2)).c = sum(C(i).sw);
    end
    if isfield(C, 'R');
        pwg(C(i).edge(1),C(i).edge(2)).R = C(i).R;
    else
        pwg(C(i).edge(1),C(i).edge(2)).R = eye(6);
    end
    if isfield(C, 'Y');
        pwg(C(i).edge(1),C(i).edge(2)).Y = C(i).Y;
    end
    if isfield(C, 'P');
        pwg(C(i).edge(1),C(i).edge(2)).P = C(i).P;
    end
    G(C(i).edge(1),C(i).edge(2)) = 1;
end