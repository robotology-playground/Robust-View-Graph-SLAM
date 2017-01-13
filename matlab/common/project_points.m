function project_points(C,pose,kpts,options)
for k=1:length(C)
        i=C(k).edge(1); % get node numbers
        j=C(k).edge(2);
        if ~isempty(C(k).xf)
            xi=kpts{i}(C(k).matches(1,:),:); % get 2d measurements
            xi=calibrate_image_points(xi,options,i); % calibrate 2d measurements
	       range=1./C(k).xf;
               X = get_scan_from_range(xi', range');
               X = transform_to_global_w(X', pose(6*(i-1)+(1:6),1))';
               plot3(X(1,:),X(2,:),X(3,:),'+','markersize',1,'color',[rand rand rand]);
	       axis equal; grid on; hold on;
	       pause;
        end
    end

%
%
function x=calibrate_image_points(x,options,k)
[K,kc]=get_intrinsics(options,k); % remove lense distortion from the key-points
x(:,1)=(x(:,1)-K(1,3))/K(1,1);
x(:,2)=(x(:,2)-K(2,3))/K(2,2);
x(:,1:2)=remove_lens_distortion(x(:,1:2),kc);

%
%
function [p, dp] = get_scan_from_range(p, r)
if size(p,1) < 3; p = pextend(p); end;
rim = sqrt(p(1,:).*p(1,:) + p(2,:).*p(2,:) + p(3,:).*p(3,:));
d = r./rim;
p(1,:) = d.*p(1,:);
p(2,:) = d.*p(2,:);
p(3,:) = d.*p(3,:);
if nargout>1
    dd = -r.*d;
    dp = p;
    dp(1,:) = dd.*p(1,:);
    dp(2,:) = dd.*p(2,:);
    dp(3,:) = dd.*p(3,:);
end
