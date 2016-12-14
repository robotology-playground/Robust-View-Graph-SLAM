function r=homography_residual(H)

% calculates the residuals for homography computation

u=H(1)*X1(1,inliers)+H(4)*X1(2,inliers)+H(7);
v=H(2)*X1(1,inliers)+H(5)*X1(2,inliers)+H(8);
d=H(3)*X1(1,inliers)+H(6)*X1(2,inliers)+1;
du=X2(1,inliers)-u./d;
dv=X2(2,inliers)-v./d;
r=sum(du.*du+dv.*dv) ;


end