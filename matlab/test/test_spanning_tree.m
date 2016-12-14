clc;
clear all;
G = round(10*rand(4,4)); 
G = G - diag(diag(G));
G = G+G';

% the algorithm
maxinl = 0;
threshold = length(G);
for jj = 1:5
    
    sw = zeros(length(G),1);
    i = ceil(rand*length(sw));
    sw(i) = 1;
    campairs = [];
    
    % Linear elimination along the edges of the tree
    while sum(sw) < length(sw);
        
        % choose the camera with most matches
        [i,j,W] = find(sw*(1-sw)'.*G);
        if isempty(W); break; end
        
        % maximum spanning tree
        if 1
            ind = find(max(W) == W);
            ind = ind(1);
            i = i(ind);
            j = j(ind);
            
        % probability weighted spanning tree
        else
            prob = rand;
            %ind = find(cumsum(W.^5)./sum(W.^5) >= prob);
            ind = find(cumsum(W)./sum(W) >= prob);
            ind = ind(1);
            i = i(ind);
            j = j(ind);
        end
        
        campairs = [campairs [i; j; W(ind)]];
        sw(j) = 1;
        
    end
    
    campairs 
    
%     res = zeros(size(G));
%     weightres = zeros(size(G));
%     inl = zeros(size(G));
%     
%     % check for inlier rotations using all the constraints
%     for i = 1:length(G);
%         for j = i+1:length(G);
%             res(i,j) = 1; % some sort of error measure that describes how close the estimated value is to the contraint
%             weightres(i,j) = G(i,j)*res(i,j);
%             if res(i,j) < threshold; inl(i,j) = 1; end
%         end
%     end
%     
%     if sum(sum(inl)) > sum(sum(maxinl));
%         maxinl = inl;
%         maxres = res;
%         maxweightres = weightres;
%         maxcampairs = campairs;
%     end

end


