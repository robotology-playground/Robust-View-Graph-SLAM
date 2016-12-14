function xs = pose_generate_sequential(C)
% Simple spanning tree: pose sequence order

edge = vertcat(C.edge);
N = max(max(edge));
xs = [0; 0; 0; 0; 0; 0];
i = 1;
done(i) = 1;
while sum(done) < N
    for k = 1:length(C)
        if C(k).edge(1) == i && C(k).edge(2) == i+1
            xs(:,i+1) = transform_to_global_w(C(k).z, xs(:,i));
            done(i+1) = 1;
            i = i+1;
        end
        sum(done)
    end
end

% xs = [0;0;0;0;0;0];
% i = 1;
% for k = 1:length(C)
%     if C(k).edge(1) == i && C(k).edge(2) == i+1
%        xs(:,i+1) = transform_to_global_w(C(k).z, xs(:,i));
%        i = i+1;
%     end
% end
