function [C, idx] = assign_constraint_weight(C)

config_rswitch_v2;

% look at disconnected nodes
edges = vertcat(C.edge);
t = zeros(max(edges(:)), 1);
t(unique(edges(:))) = 1;
idx = find(t==0)
% remove disconnected nodes
for i=1:length(C)
    C(i).edge(1) = C(i).edge(1) - sum(C(i).edge(1)>idx);
    C(i).edge(2) = C(i).edge(2) - sum(C(i).edge(2)>idx);
end

% convert from constraints to graph
C = apply_nominal_R(C);
[C, G] = C_to_pwg(C);

% look at the number of constraints per node
thickness = sum(G, 2);
if 0; figure; plot(thickness); end;

switch EDGE_WEIGHT;
    case 1
        % loop statistics
        disp('    ');
        disp('Gathering loop statistics.');
        C = constraint_loop_statistics(C);
    case 2
        % number of tracks
        disp('    ');
        disp('Weights by number of tracks.');
        C = constraint_tracks(C);
    case 3
        % trace of information
        disp('    ');
        disp('Weights by trace of information.');
        C = constraint_information(C);
    case 4
        % trace of covariance
        disp('    ');
        disp('Weights by trace of covariance.');
        C = constraint_covariance(C);
end

% convert from graph to constraints
C = convert_C_to_CG(C);
%CG = pwg_to_C(C);
%
%
function C = constraint_tracks(C)
for i = 1:size(C,1)
    for j = (i+1):size(C,2)
        if isempty(C(i,j).z), continue, end
        C(i,j).yes = 0;
        C(i,j).no = 0;
    end
end
for i = 1:size(C,1)
    for j = (i+1):size(C,2)
        if isempty(C(i,j).z), continue, end
        if isfield('C','c'); C(i,j).yes = C(i,j).c; end
        if isfield('C','w'); C(i,j).yes = C(i,j).w; end
    end
end
%
%
function C = constraint_information(C)
for i = 1:size(C,1)
    for j = (i+1):size(C,2)
        if isempty(C(i,j).z), continue, end
        C(i,j).yes = 0;
        C(i,j).no = 0;
    end
end
for i = 1:size(C,1)
    for j = (i+1):size(C,2)
        if isempty(C(i,j).z), continue, end
        C(i,j).yes = trace(C(i,j).Y);
    end
end
%
%
function C = constraint_covariance(C)
for i = 1:size(C,1)
    for j = (i+1):size(C,2)
        if isempty(C(i,j).z), continue, end
        C(i,j).yes = 0;
        C(i,j).no = 0;
    end
end
for i = 1:size(C,1)
    for j = (i+1):size(C,2)
        if isempty(C(i,j).z), continue, end
        C(i,j).yes = -trace(C(i,j).P);
    end
end