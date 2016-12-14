% Returns indexes of 3D line segment pairs {[X0 Y0],[X Y]} violating
% ordering constraint wrt camera centers C1, C2.
% The result is invariant to transforming its parameters by any homography.

function n = ordering_constraint_pairwise_3d(C1,C2,X0,Y0,X,Y)

B = vgg_contreps(C1*C2'-C2*C1');  % dual of line joining camera centres

% Make det([C1 C2 X Y]) positive for all segments,
% ie, swap segments that rotate differently wrt baseline than [X0 Y0].
if X0'*B*Y0 < 0  % note: X0'*B*Y0 = det([C1 C2 X0 Y0])
  Z=X0; X0=Y0; Y0=Z;
end
n = sum(X.*(B*Y),1) < 0;
Z=X(:,n); X(:,n)=Y(:,n); Y(:,n)=Z;

% Find [Xn Yn] segments that are intersected by epipolar beam associated with [X0 Y0]
n = find( (X0'*B*Y>0) & (Y0'*B*X<0) ); % segments intersected by the pencil of epipolar planes spanned by [X0 Y0]
Xn = X(:,n);
Yn = Y(:,n);

% Find X1 intersections of lines [X Y] with epipolar plane [C1 C2 X0]
E = X0'*B; % epipolar plane passing thru X0
EX = E*Xn;
EY = E*Yn;
X1 = EX([1 1 1 1],:).*Yn - EY([1 1 1 1],:).*Xn; % intersections of E with line segments [X Y]

% Find X2 intersections of lines [X Y] with epipolar plane [C1 C2 X0]
E = Y0'*B; % epipolar plane passing thru Y0
EX = E*Xn;
EY = E*Yn;
Y1 = EX([1 1 1 1],:).*Yn - EY([1 1 1 1],:).*Xn; % intersections of E with line segments [X Y]

% Ordering constraint is expressed as follows:
% Simplices [C1 X0 Z U] and [C2 X0 Z U] must have the same orientation,
% for any auxilliary point U off the plane [C1 C2 X0]. We chose U=Y0.
CXY1 = vgg_wedge([C1 X0 Y0]);
CXY2 = vgg_wedge([C2 X0 Y0]);
n = n( (CXY1*X1).*(CXY2*X1)<0 & (CXY1*Y1).*(CXY2*Y1)<0 );

return