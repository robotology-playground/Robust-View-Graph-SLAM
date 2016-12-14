function [t, nk] = minimum_spanning_tree(i, j, w)
% Compute minimum-weight spanning tree according to Kruskal's algorithm.
% Adapted from minspan.m by Michael G. Kay, which was obtained from Kevin
% Murphy's graph toolbox.

N = length(w);
assert(N > 0, 'Graph cannot be empty')
assert(length(i) == N && length(j) == N)

% Order for minimum weight
[~, iw] = sort(w);
i = i(iw);
j = j(iw);

% Compute tree
M = max([max(i) max(j)]);
v = 1:M;        % arc labels
t = zeros(1, N);% arcs in spanning tree
nt = 0;         % number of arcs in spanning tree
for k = 1:N     % k is current arc
    vi = v(i(k));
    vj = v(j(k));
    if vi ~= vj
        v(v==vj) = vi; % FIXME: this is inefficient
        t(k) = 1;
        nt = nt + 1;
        if nt == M-1, break, end
    end
end
iwr = inverse_mapping(iw);
t = t(iwr); 

% Check whether graph is connected, otherwise we have a spanning forest
if nargout == 2
    c = unique(v([i; j]));      % unique labels of arc vertices 
    %c = unique(v);              
    nk = length(c);
    if any(t)==0, nk = 0; end   % self-loop is not a component
end

% Note: c = unique(v([i; j])) detects the unique labels of only vertices
% that appear in {i,j}. This might not include all vertices in 1:M. The
% alternative, c = unique(v), checks the labels of all vertices in 1:M. 

%
%

function imap = inverse_mapping(map)
imap = zeros(size(map));
imap(map) = 1:length(map);
