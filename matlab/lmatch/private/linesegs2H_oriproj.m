% u1 ... segment's end points in image 1
% u2 ... segment's end points in image 2
% u1, u2 are two *exactly corresponding* point pairs
function H = linesegs2H_oriproj(u1,u2,P1,P2)

F = vgg_F_from_P(P1,P2);
e2 = P2*vgg_wedge(P1);

F(:) = normx(F(:));
e2 = normx(e2);

% Family of homographies consistent with P1,P2 is H(m) = H0 + e2*m, 1-by-3 vector m is the free parameter.
H0 = vgg_contreps(e2)*F;

% Find s such that it must be m*hom(u1) = s
% From s, 1-parameter family of m is m(t) = m0 + t*l1.
for n = 1:2
  a(:,n) = vgg_contreps([u2(:,n);1])*H0*[u1(:,n);1];
  b(:,n) = vgg_contreps([u2(:,n);1])*e2;
  s(n) = -b(:,n)\a(:,n);
end
m0 = s/hom(u1);
l1 = vgg_wedge([u1; 1 1]);

% Update H0 so that family of homographies consistent with P1,P2 and
% sending u1 to u2 is H(t) = H0 + t*e2*l1.
H0 = H0 + e2*m0;

% The parameter t is found such that it makes H preserve area of
% a surface element in the middle of the line segment.
% That means, it sets Jacobian of the mapping by u := nhom(H*hom(u)) to 1.
%
% This Jacobian turns out to be linear function of t, as long as the
% area element is on the line segment, l1*hom(u)==0 because e2*l1 vanishes.

% compute Jacobian of mapping u' = nhom(H*hom(u)) in point u=u1(:,1) is J(t) = det(M0 + t*Mt)
M = [eye(2) -u2(:,1)]/(H0(3,:)*[u1(:,1);1]);
M0 = M*H0(:,1:2);
Mt = M*e2*l1(1:2);

% Rewrite J(t) using det(A+B)=det(A)+det(B)+det([A(:,1) B(:,2)])+det([B(:,1) A(:,2)]).
% Since det(Mt)==0, it is J(t) = det(M0) + a*t.
a = det([M0(:,1) Mt(:,2)])+det([Mt(:,1) M0(:,2)]);
t = (1-det(M0))/a;

H = H0 + t*e2*l1;

return
