function [A, B] = test_linear_system(P1, P2, U, u1, u2)

% to test, run the lines:
% numpts=4;  % number of 3d points
% [A, B] = test_linear_system(rand(3, 4), rand(3, 4), rand(4, numpts), rand(3, numpts), rand(3, numpts)); 
% spy(A);

% rotation specific vectors
Ba = [0 1 0; -1 0 0; 0 0 0];
Bb = [0 0 1; 0 0 0; -1 0 0];
Bc = [0 0 0; 0 0 1; 0 -1 0];

% allocate vectors for the system matrices
resnum = 0;
resnr = size(u1, 2);
row = zeros(2*resnr*(6+3), 1);
col = zeros(2*resnr*(6+3), 1);
data = zeros(2*resnr*(6+3), 1);
B = zeros(2*resnr, 1);
lastentry = 0;
lastentryB = 0;

% prepair data points
U = pflat(U);
U = U(1:3,:);
P = {P1, P2};
u = [u1; u2];
numpts = size(u1, 2);
vis = 1:numpts;

for i = 1:2;
    
    % initial camera parameters
    R0 = P{i}(:,1:3);
    t0 = P{i}(:,4);
    
    % calculate derivatives for both residuals in all images
    % U1,U2,U3 - 3d point parameters
    % a,b,c - rotation parameters for camera i
    % t1,t2,t3 - translation parameters for camera i.
    % 3d point derivatives
    
    % common factors
    num1 = R0(1,:)*U + t0(1);
    num2 = R0(2,:)*U + t0(2);
    denum1 = R0(3,:)*U + t0(3);
    denum2 = denum1.^2;
    
    % the derivatives
    dudt1 = 1./denum1;
    dvdt1 = zeros(size(dudt1));
    dudt2 = zeros(size(dudt1));
    dvdt2 = 1./denum1;
    dudt3 = -num1./denum2;
    dvdt3 = -num2./denum2;
    
    duda = (Ba(1,:).*R0(1,:)*U)./denum1 - num1./denum2.*(Ba(3,:).*R0(3,:)*U);
    dvda = (Ba(2,:).*R0(2,:)*U)./denum1 - num2./denum2.*(Ba(3,:).*R0(3,:)*U);
    dudb = (Bb(1,:).*R0(1,:)*U)./denum1 - num1./denum2.*(Bb(3,:).*R0(3,:)*U);
    dvdb = (Bb(2,:).*R0(2,:)*U)./denum1 - num2./denum2.*(Bb(3,:).*R0(3,:)*U);
    dudc = (Bc(1,:).*R0(1,:)*U)./denum1 - num1./denum2.*(Bc(3,:).*R0(3,:)*U);
    dvdc = (Bc(2,:).*R0(2,:)*U)./denum1 - num2./denum2.*(Bc(3,:).*R0(3,:)*U);
    
    dudx = R0(1,1)./denum1 - num1./denum2.*R0(3,1);
    dvdx = R0(2,1)./denum1 - num2./denum2.*R0(3,1);
    dudy = R0(1,2)./denum1 - num1./denum2.*R0(3,2);
    dvdy = R0(2,2)./denum1 - num2./denum2.*R0(3,2);
    dudz = R0(1,3)./denum1 - num1./denum2.*R0(3,3);
    dvdz = R0(2,3)./denum1 - num2./denum2.*R0(3,3);
    
    if 0
        aa = 2*rand(1,1)-1;
        bb = 2*rand(1,1)-1;
        cc = 2*rand(1,1)-1;
        UU1 = 2*rand(1,1)-1;
        UU2 = 2*rand(1,1)-1;
        UU3 = 2*rand(1,1)-1;
        tt1 = 2*rand(1,1)-1;
        tt2 = 2*rand(1,1)-1;
        tt3 = 2*rand(1,1)-1;
        
        f = [];
        kk = ceil(rand*length(vis));
        kk = vis(kk);
        f01 = (-u(2*i-1,kk)+(R0(1,:)*U(:,kk)+ t0(1))./(R0(3,:)*U(:,kk)+t0(3)));
        f02 = (-u(2*i,kk)+(R0(2,:)*U(:,kk)+ t0(2))./(R0(3,:)*U(:,kk)+ t0(3)));
        
        f2 = [];
        for s = -0.1:0.001:0.1
            R = expm(s*(Ba*aa+Bb*bb+Bc*cc))*R0;
            t = t0+s*[tt1;tt2;tt3];
            UU = U+s*repmat([UU1;UU2;UU3], [1 size(U,2)]);
            f = [f; (u(2*i-1,kk)-(R(1,:)*UU(:,kk)+t(1))./(R(3,:)*UU(:,kk)+ t(3)))^2 + ...
                (u(2*i,kk)-(R(2,:)*UU(:,kk)+t(2))./(R(3,:)*UU(:,kk)+ t(3)))^2];
            f2 = [f2; (f01 + s*(duda(kk)*aa + dudb(kk)*bb + dudc(kk)*cc + ...
                dudx(kk)*UU1 + dudy(kk)*UU2 + dudz(kk)*UU3 + ...
                dudt1(kk)*tt1 + dudt2(kk)*tt2 + dudt3(kk)*tt3)).^2 + ...
                (f02 + s*(dvda(kk)*aa + dvdb(kk)*bb + dvdc(kk)*cc + ...
                dvdx(kk)*UU1 + dvdy(kk)*UU2 + dvdz(kk)*UU3 + ...
                dvdt1(kk)*tt1 + dvdt2(kk)*tt2 + dvdt3(kk)*tt3)).^2];
        end
        %figure(1);plot(-0.1:0.001:0.1,[f f2(end:-1:1)]);
        figure(4); plot(-0.1:0.001:0.1,f);
        figure(5); plot(-0.1:0.001:0.1,f2);
    end
    
    % collect first residuals
    % 3d point parameters:
    % U1-coeffs
    row( lastentry + (1:numpts) ) = resnum + (1:2:2*numpts)';
    col( lastentry + (1:numpts) ) = ( ((1:numpts)-1)*3 + 1 )';
    data( lastentry + (1:numpts) ) = dudx';
    lastentry = lastentry + numpts;
    % U2-coeffs
    row( lastentry + (1:numpts) ) = resnum + (1:2:2*numpts)';
    col( lastentry + (1:numpts) ) = ( ((1:numpts)-1)*3 + 2 )';
    data( lastentry + (1:numpts) ) = dudy';
    lastentry = lastentry + numpts;
    % U3-coeffs
    row( lastentry + (1:numpts) ) = resnum + (1:2:2*numpts)';
    col( lastentry + (1:numpts) ) = ( (1:numpts)*3 )';
    data( lastentry + (1:numpts) ) = dudz';
    lastentry = lastentry + numpts;
    % camera parameters:
    % a-coeffs
    row( lastentry + (1:numpts) ) = resnum + (1:2:2*numpts)';
    col( lastentry + (1:numpts) ) = (3*numpts + (i-1)*6 + 1)*ones(numpts, 1);
    data( lastentry + (1:numpts) ) = duda';
    lastentry = lastentry + numpts;
    % b-coeffs
    row( lastentry + (1:numpts) ) = resnum + (1:2:2*numpts)';
    col( lastentry + (1:numpts) ) = (3*numpts + (i-1)*6 + 2)*ones(numpts, 1);
    data( lastentry + (1:numpts) ) = dudb';
    lastentry = lastentry + numpts;
    % c-coeffs
    row( lastentry + (1:numpts) ) = resnum + (1:2:2*numpts)';
    col( lastentry + (1:numpts) ) = (3*numpts + (i-1)*6 + 3)*ones(numpts, 1);
    data( lastentry + (1:numpts) ) = dudc';
    lastentry = lastentry + numpts;
    % t_1-coeffs
    row( lastentry + (1:numpts) ) = resnum+ (1:2:2*numpts)';
    col( lastentry + (1:numpts) ) = (3*numpts + (i-1)*6 + 4)*ones(numpts, 1);
    data( lastentry + (1:numpts) ) = dudt1';
    lastentry = lastentry+numpts;
    % t_2-coeffs
    row( lastentry + (1:numpts) ) = resnum+ (1:2:2*numpts)';
    col( lastentry + (1:numpts) ) = (3*numpts + (i-1)*6 + 5)*ones(numpts, 1);
    data( lastentry + (1:numpts) ) = dudt2';
    lastentry = lastentry + numpts;
    % t_3-coeffs
    row( lastentry + (1:numpts) ) = resnum+ (1:2:2*numpts)';
    col( lastentry + (1:numpts) ) = (3*numpts + (i-1)*6 + 6)*ones(numpts, 1);
    data( lastentry + (1:numpts) ) = dudt3';
    lastentry = lastentry + numpts;
    
    % collect second parameters:
    % 3d point parameters:
    % U1-coeffs
    row( lastentry + (1:numpts) ) = resnum + (2:2:2*numpts)';
    col( lastentry + (1:numpts) ) = ( ((1:numpts)-1)*3 + 1 )';
    data( lastentry + (1:numpts) ) = dvdx';
    lastentry = lastentry + numpts;
    % U2-coeffs
    row( lastentry + (1:numpts) ) = resnum + (2:2:2*numpts)';
    col( lastentry + (1:numpts) ) = ( ((1:numpts)-1)*3 + 2 )';
    data( lastentry + (1:numpts) ) = dvdy';
    lastentry = lastentry+numpts;
    % U3-coeffs
    row( lastentry + (1:numpts) ) = resnum + (2:2:2*numpts)';
    col( lastentry + (1:numpts) ) = ( ((1:numpts)-1)*3 + 3 )';
    data( lastentry + (1:numpts) ) = dvdz';
    lastentry = lastentry + numpts;
    % camera parameters:
    % a-coeffs
    row( lastentry + (1:numpts) ) = resnum + (2:2:2*numpts)';
    col( lastentry + (1:numpts) ) = (3*numpts + (i-1)*6 + 1)*ones(numpts, 1);
    data( lastentry + (1:numpts) ) = dvda';
    lastentry = lastentry + numpts;
    % b-coeffs
    row( lastentry + (1:numpts) ) = resnum + (2:2:2*numpts)';
    col( lastentry + (1:numpts) ) = (3*numpts + (i-1)*6 + 2)*ones(numpts, 1);
    data( lastentry + (1:numpts) ) = dvdb';
    lastentry = lastentry + numpts;
    % c-coeffs
    row( lastentry + (1:numpts) ) = resnum + (2:2:2*numpts)';
    col( lastentry + (1:numpts) ) = (3*numpts + (i-1)*6 + 3)*ones(numpts, 1);
    data( lastentry + (1:numpts) ) = dvdc';
    lastentry = lastentry + numpts;
    % t_1-coeffs
    row( lastentry + (1:numpts) ) = resnum + (2:2:2*numpts)';
    col( lastentry + (1:numpts) ) = (3*numpts + (i-1)*6 + 4)*ones(numpts, 1);
    data( lastentry + (1:numpts) ) = dvdt1';
    lastentry = lastentry + numpts;
    % t_2-koeff
    row( lastentry + (1:numpts) ) = resnum + (2:2:2*numpts)';
    col( lastentry + (1:numpts) ) = (3*numpts + (i-1)*6 + 5)*ones(numpts, 1);
    data( lastentry + (1:numpts) ) = dvdt2';
    lastentry = lastentry + numpts;
    % t_3-coeffs
    row( lastentry + (1:numpts) ) = resnum + (2:2:2*numpts)';
    col( lastentry + (1:numpts) ) = (3*numpts + (i-1)*6 + 6)*ones(numpts, 1);
    data( lastentry + (1:numpts) ) = dvdt3';
    lastentry = lastentry + numpts;
    
    % increment residuals counter
    resnum = resnum + 2*numpts;
    
    % constant terms
    btmp = zeros(2*numpts, 1);
    % first residuals
    btmp(1:2:end) = (P{i}(1,:)*pextend(U))./(P{i}(3,:)*pextend(U))-u(2*i-1,:);
    % second residuals
    btmp(2:2:end) = (P{i}(2,:)*pextend(U))./(P{i}(3,:)*pextend(U))-u(2*i,:);
    % B = [B; btmp];
    B( lastentryB + (1:length(btmp)) ) = btmp;
    lastentryB = lastentryB + length(btmp);
    
end

% fill in all parameters
A = sparse(row,col,data);

% Lock the coordinate system by setting the translation and rotation terms
% of the first camera to zero
%A = A(:, [1:3*numpts 3*numpts+7:end]);
%A = A(:, 2:end);  % for some reason, the x coordinate of the first point was removed