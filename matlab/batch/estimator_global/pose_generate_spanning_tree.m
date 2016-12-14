function [x, t] = pose_generate_spanning_tree(C, t)
config_rswitch_v2;
if nargin == 1
    edges = vertcat(C.edge);
    w = -[C.w]; % if using number of matches or inliers
    [t, nk] = minimum_spanning_tree(edges(:,1), edges(:,2), w);
    assert(nk == 1,'Tree must be connected')
end

Ct = C(t==1);
N = length(Ct);
x = zeros(6,N+1);%x(1,1)=1;

doneC = false(1, N);
donepose = [true false(1, N)]; 
while 1
    alldone = 1;
    for i = 1:N
        if doneC(i) == 1,continue,end
        Ci = Ct(i);
        i1 = Ci.edge(1);        
        i2 = Ci.edge(2);        
        
        if donepose(i2) == 1     % i2 only is known, swap (i1,i2)
            Ci.z = transform_to_relative_w([0;0;0;0;0;0], Ci.z);
            tmp = i1; i1 = i2; i2 = tmp;            
        elseif donepose(i1) == 0 % i1 and i2 unknown
            alldone = 0;
            continue
        end
        
        assert(donepose(i1) == 1 && donepose(i2) == 0) % i1 only is known
        x(:, i2) = transform_to_global_w(Ci.z, x(:, i1));
        %x(i2,1:4)=forward_rotation(Ci.q,x(i1,1:4));
        %x(i2,5:7)=forward_position(Ci.t,x(i1,5:7),x(i1,1:4));
        %p=[Ci.q Ci.t]';
        %b=x(i1,:)';
        %p=transform_to_global(p,b);
        %x(i2,:)=p(:,1)';
        donepose(i2) = 1;
        doneC(i) = 1;
    end
    if alldone, break, end
end

if 0
    plot3(x(1,:), x(3,:), x(2,:),'.:'), hold on
    for i = 1:N
        plot3(x(1, Ct(i).edge), x(3, Ct(i).edge), x(2, Ct(i).edge), 'g')
    end
    axis equal
end