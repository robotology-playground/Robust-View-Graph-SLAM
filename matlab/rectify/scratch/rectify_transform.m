function [H1,H2,newe,G2,R2] = rect_transform(P1, P2, e1, e2, x1, x2)

%
%   Input arguments are:
%   - P1 and P2 should both be 3-by-4 projection matrices, with P1 = [I 0].
%   - e1 and e2 should both be 3-by-1 vectors containing the two epipoles, with e1
%     being the epipole in image 1 (the projection of the optical centre of the
%     second camera onto the first image), and e2 being the epipole in image 2
%     (projection of the optical centre of the first camera onto the second image).
%   - x1 and x2 should both be 3-by-n matrices with each column of the matrix being
%     an image point in homogeneous coordinates and n being is the number of matching
%     points.
%
%   Output arguments are:
%   - H1 and H2 would both be 3-by-3 rectification matrices.
%   - G2 and R2 are all optional output arguments.  If specified, they would
%     all be 3-by-3 matrices, satisfying the condition H2 = G2*R2.
%
% rect_transform.m : computes the rectification transformation matrices
% section 11.12 Hartley

% check input argument P1 -- it must be a [I 0] matrix
if max(max(abs(P1-[eye(3) zeros(3,1)]))) ~= 0
    error('rectification_transf: matrix P1 must be of the form [I 0]');
end

% compute rectification matrix H2 for the second image
if nargout >= 4 & nargout <= 5
    [H2,newe,G2,R2] = rectify_H(e2);
else
    [H2,newe] = rectify_H(e2);
end

H0 = H2*P2(:,1:3);
x2hat = pflat(H2*x2);
x1hat = pflat(H0*x1);
B = [
    sum( (ones(3,1)*x1hat(1,:)) .* x1hat, 2 )';
    sum( (ones(3,1)*x1hat(2,:)) .* x1hat, 2 )';
    sum( (ones(3,1)*x1hat(3,:)) .* x1hat, 2 )'
    ];

if sum(abs(newe-[1;0;0])) == 0
    b = [
        sum( x1hat(1,:) .* (x2hat(1,:)-x1hat(1,:)) );
        sum( x1hat(2,:) .* (x2hat(1,:)-x1hat(1,:)) );
        sum( x1hat(3,:) .* (x2hat(1,:)-x1hat(1,:)) )
        ];
else
    b = [
        sum( x1hat(1,:) .* (x2hat(2,:)-x1hat(2,:)) );
        sum( x1hat(2,:) .* (x2hat(2,:)-x1hat(2,:)) );
        sum( x1hat(3,:) .* (x2hat(2,:)-x1hat(2,:)) )
        ];
end

abc = B \ b;

if sum(abs(newe-[1;0;0])) == 0
    A = [[1 0 0]+abc'; 0 1 0; 0 0 1];
else
    A = [1 0 0; [0 1 0]+abc'; 0 0 1];
end
H1 = A*H0;

return


