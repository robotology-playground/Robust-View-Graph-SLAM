% Finds incompatible segments in 2 images due to ordering constraint.
%   i ... indices of compatible segments
function i = ordering_constraint_pairwise(li0,li,l,P,I)
% The result is invariant to transforming P by any 3-D homography.

if isempty(li)
  i = [];
  return
end

F = vgg_F_from_P(P{:});  F = F/sqrt(sum(F(:).^2));
e1 = normx(P{1}*vgg_wedge(P{2}));
e2 = normx(P{2}*vgg_wedge(P{1}));

% Inititialize aux. variables:
%   x11,x12,l1 ... start point / end point / straight line of segment 1 in image 1
%   x21,x22,l2 ... start point / end point / straight line of segment 1 in image 2
%   y11,y12,m1 ... start points / end points / straight lines of segments 2 in image 1
%   y21,y22,m2 ... start points / end points / straight lines of segments 2 in image 2
s1 = l{1}(li0(1));   x11 = [s1.u;1];              x12 = [s1.v;1];              l1 = s1.l;
s2 = l{2}(li0(2));   x21 = [s2.u;1];              x22 = [s2.v;1];              l2 = s2.l;
t1 = l{1}(li(1,:));  y11 = [t1.u]; y11(3,:) = 1;  y12 = [t1.v]; y12(3,:) = 1;  m1 = vertcat(t1.l);
t2 = l{2}(li(2,:));  y21 = [t2.u]; y21(3,:) = 1;  y22 = [t2.v]; y22(3,:) = 1;  m2 = vertcat(t2.l);

% Swap directions of all segments so that - 
%   -  l1*e1 > 0, m1*e1 > 0  in image 1
%   -  l2*e2 < 0, m2*e2 < 0  in image 2
[x11,x12,l1] = swap_direction( +e1, x11,x12,l1 );
[x21,x22,l2] = swap_direction( -e2, x21,x22,l2 );
[y11,y12,m1] = swap_direction( +e1, y11,y12,m1 );
[y21,y22,m2] = swap_direction( -e2, y21,y22,m2 );
%figure(1); clf;  plotuvo(nhom(y11),nhom(y12),'y'); axis image; hold on; figs full
%figure(2); clf;  plotuvo(nhom(y21),nhom(y22),'y'); axis image; hold on; figs full

i = [];

if isempty( beam_search( -(F *[x11 -x12])', x22,x21 ) ) % Check if the 1st segments interact at all
  return
end

% Clip 1st segment by epi. beam
[x22,x21] = beam_clip( -(F *[x11 -x12])', x22,x21,-l2 );
if isempty(x22) | isempty(x21)
  return
end
[x12,x11] = beam_clip( +(F'*[x21 -x22])', x12,x11,-l1 );
if isempty(x12) | isempty(x11)
  return
end
% a=images(I);
% axes(a(1)); plotuvo(nhom(x11),nhom(x12),'r');  plotl((F'*[x21 x22])','b');
% axes(a(2)); plotuvo(nhom(x21),nhom(x22),'r');  plotl((F *[x11 x12])','b');

% Find 2nd segments interacting with 1st segment
b1 = (F'*[x21 -x22])';  i = beam_search( +b1, y12,y11 );
b2 = (F *[x11 -x12])';  i = i( beam_search( -b2, y22(:,i),y21(:,i) ) );
y11 = y11(:,i);  y12 = y12(:,i);  m1 = m1(i,:);
y21 = y21(:,i);  y22 = y22(:,i);  m2 = m2(i,:);
% axes(a(1)); plotuvo(nhom(y11),nhom(y12),'g');
% axes(a(2)); plotuvo(nhom(y21),nhom(y22),'g');

% Clip 2nd segments with 1st segments
[y12,y11] = beam_clip( +b1, y12,y11,-m1 );
if isempty(y12) | isempty(y11)
  return
end
[y22,y21] = beam_clip( -b2, y22,y21,-m2 );
if isempty(y22) | isempty(y21)
  return
end
% axes(a(1)); plotuvo(nhom(y11),nhom(y12),'y');
% axes(a(2)); plotuvo(nhom(y21),nhom(y22),'y');

% Clip 1st segment by 2nd segments
ce1 = vgg_contreps(e1);
ce2 = vgg_contreps(e2);
p = vgg_contreps(l1)*ce1;  x11 = p*y11;  x12 = p*y12;
p = vgg_contreps(l2)*ce2;  x21 = p*y21;  x22 = p*y22;
%Test: all these numbers must be zero:
%norm( wedge(normx(x11),normx(y11))*normx(e1) )
%norm( wedge(normx(x12),normx(y12))*normx(e1) )
%norm( wedge(normx(x21),normx(y21))*normx(e2) )
%norm( wedge(normx(x22),normx(y22))*normx(e2) )

i = i( ( sum(cros(x11,y11).*(ce1*x11)) .*...
         sum(cros(x21,y21).*(ce2*x21)) > 0 ) |...
       ( sum(cros(x12,y12).*(ce1*x12)) .*...
         sum(cros(x22,y22).*(ce2*x22)) > 0 ) );

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Swaps direction of segments (x1,x2,l) such that l*e>0.
function [x1,x2,l] = swap_direction(e,x1,x2,l)
n = l*e < 0;
l(n,:) = -l(n,:);
aux = x1(:,n);  x1(:,n) = x2(:,n);  x2(:,n) = aux;
return


% n = beam_search(b,x1,x2)  Finds indices of segments (x1,x2) that intersect beam b.
%
% Beam b is space between 2 lines, b = [l1;l2]. Ie, set of points x such that b*x>0.
% x1,x2 (in homog. coords.) are assumed :-
%   - to be such that wedge(x1,x2)*wedge(b) > 0,
% n satisfies: all( b(1,:)*x1(:,n)>0 & b(2,:)*x2(:,n)>0 ).
function n = beam_search(b,x1,x2)

n = 1:size(x1,2);

n( b(1,:)*x1(:,n) < 0 ) = []; % remove segments with end point above first epipolar 
if isempty(n)
  return
end

n( b(2,:)*x2(:,n) < 0 ) = []; % remove segments with start point under second epipolar
return


% [x1,x2] = beam_clip(b,x1,x2 [,l])  Clips segments (x1,x2) by beam b.
%
% Beam b is space between 2 lines, b = [l1;l2]. Ie, set of points x such that b*x>0.
% x1,x2 (in homog. coords.) are assumed :-
%   - to be such that wedge(x1,x2)*wedge(b) > 0,
%   - to have nonempty intersection with beam b, ie, all( b(1,:)*x1>0 & b(2,:)*x2>0 ).
% This is preserved at the output.
% l is auxiliary to speed things up, l(i,:)=wedge([x1(:,i) x2(:,i)]).
function [x1,x2] = beam_clip(b, x1,x2,l)

n = b(2,:)*x1 < 0;  
x1(:,n) = vgg_contreps(+b(2,:)) * l(n,:)';  % clip startpoint

n = b(1,:)*x2 < 0;
x2(:,n) = vgg_contreps(-b(1,:)) * l(n,:)';  % clip endpoint

return

