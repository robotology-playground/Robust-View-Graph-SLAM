function [x,P] = recover_moments(y,Y,idx,ncams)
%[x,P] = recover_moments(y,Y)
%
% Computes moment-form estimate {x,P} given an infomation-from estimate {y,Y}.
%
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2014

if ~issparse(Y)
    Y = sparse(Y);
end
x = y; % memory initialisation
P = Y; % using two variables for change
i = amd(Y); % compute reordering

warning off MATLAB:singularMatrix
warning off MATLAB:nearlySingularMatrix

try L = lchol(Y(i,i));
	disp( 'Matrix PD' );
    f = L\y(i,1); % triangular solve: L\y
    x(i,1) = full(L'\f); % mean with original ordering
    if nargout > 1
        P(i,i) = spinv(Y(i,i)); % try using Takahashi's equations 
        % file located at (~/Dev/GPStuff-4.6/gp/private)
        P = full(P);
    end
catch % reconditioning again, then try using Takahashi's equations.
    beta = 1e-8;
    disp( ['Reconditioning with beta = ',num2str(beta)] );
    try L = lchol(Y(i,i) + speye(size(Y))*beta);
        f = L\y(i,1); % triangular solve: L\y
        x(i,1) = full(L'\f); % mean with original ordering
        if nargout > 1
            P = spinv(Y);
            P = full(P);
        end
    catch
        beta = 1e-6; % was 1 before
        disp( ['Reconditioning with beta = ',num2str(beta)] );
        try L = lchol(Y(i,i) + speye(size(Y))*beta);
            f = L\y(i,1); % triangular solve: L\y
            x(i,1) = full(L'\f); % mean with original ordering
            if nargout > 1
                P = spinv(Y);
                P = full(P);
            end
        catch
            if nargin>2
                disp( 'Y is not PSD, using lchol() to computer the marginal' );
                x = zeros(length(Y),1);
                try L = lchol(Y(idx,idx));
                    f = L\y(idx,1); % triangular solve: L\y
                    x(idx,1) = full(L'\f); % mean with original ordering
                    if nargout > 1
                        if size(idx, 2)==1; idx = idx'; end;
                        j = repmat(idx, length(idx), 1);
                        i = j';
                        s = full(spinv(Y(idx,idx)));
                        P = sparse(i(:),j(:),s(:),length(Y),length(Y));
                        P = full(P);
                    end
                catch
                    disp( 'Y is not PSD, using inv() to compute the marginal' );
                    %gcp;
                    %spmd; obj = Y(idx,idx)\y(idx,1); end;
                    %x(idx,1) = full(obj{1});
                    %spmd; x = full(Y\y); end;
                    %spmd; x(idx,1) = full(Y(idx,idx)\y(idx,1)); end;
                    x(idx,1)=Y(idx,idx)\y(idx,1);
                    %x = marginal_solve_blocks(y, Y, idx, ncams);
                    if nargout > 1
                        %N = size(Y,1);
                        %N = length(idx);
                        %I=distributed.speye(N); % distribute the sparse identity matrix
                        %I = speye(N);
                        %spmd; P = full(Y\I); end;
                        %spmd; obj = Y(idx,idx)\I; end;
                        %P(idx,idx) = full(obj{1});
                        if size(idx,2)==1; idx = idx'; end;
                        P=marginal_invert_blocks(Y,idx,ncams);
                        P=full(P);
                    end
                    %delete(gcp);
                end
            else
                disp( 'Y is not PSD, using inv() to compute the full covariance' );
                x(i,1) = full(Y(i,i)\y(i,1));
                if nargout > 1
                    P = inv(Y);
                    P = full(P);
                end
            end
        end
    end
end

warning on MATLAB:singularMatrix
warning on MATLAB:nearlySingularMatrix

%
%
function x = marginal_solve_blocks(y, Y, idx, ncams)
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2014
block_size = 1;
x = zeros(length(y), 1);
% update the points starting from the second block
%for k = (6*ncams+1):block_size:length(idx)
for k = (6*ncams+block_size+1):block_size:length(idx)
    i1 = 1:6*ncams;
    i2 = k:min(k+block_size-1,length(idx));
    jj = idx([i1, i2]);
    x(jj,1) = full(Y(jj,jj)\y(jj,1));
end
% update the points in the first block and the poses
%i1 = 1:6*ncams;
%i2 = 6*ncams+(1:block_size);
%jj = idx([i1, i2]);
jj = idx(i1);
x(jj,1) = full(Y(jj,jj)\y(jj,1));
%
%
function P = marginal_invert_blocks(Y, idx, ncams)
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2014
block_size = 1;
j = repmat(idx, length(idx), 1);
i = j';
s = zeros(length(idx));
%s = full(inv(Y(idx,idx)));
% update the points starting from the second block
%for k = (6*ncams+1):block_size:length(idx)
for k = (6*ncams+1+block_size):block_size:length(idx)
    i1 = 1:6*ncams;
    i2 = k:min(k+block_size-1,length(idx));
    ii = [i1, i2];
    jj = idx([i1, i2]);
%     if condest(Y(jj,jj))>1e12
%         A = Y(jj,jj)
%         condest(Y(jj,jj))
%         pause
%     end
    s(ii,ii) = full(inv(Y(jj,jj)));
end
% update the points in the first block and the poses
i1 = 1:6*ncams;
i2 = 6*ncams+(1:block_size);
ii = [i1, i2];
jj = idx([i1, i2]);
s(ii,ii) = full(inv(Y(jj,jj)));
P = sparse(i(:),j(:),s(:),length(Y),length(Y));

% try [LD,p,q] = ldlchol(Y,1e-8);
%   f = LD\y(q,1); % triangular solve: L\y
%	x(q,1) = full(LD'\f); % mean with original ordering
%	I = speye(size(Y));
%	F = full(L\I);
%	P(:,q) = full(U)\F;
% catch
%         end
%         [L,U,p] = lu(Y(i,i),'vector');
%         disp('using LU-Factorisation');
%         %x = U\(L\(y(p)));
%         y = y(i,1);
%         f = L\y(p); % triangular solve: L\y
%         x(i,1) = full(U\f); % mean with original ordering
%         if nargout > 1
%             F = L\speye(size(Y));
%             P(:,p) = (U\F);
%         end

%         Y0 = diag(diag(Y));
%         Y( Y > min(diag(Y)) ) = 0;
%         Y = Y + Y0;
%         Y = (Y + Y')/2;
%         i = amd(Y); % compute reordering
%         L = lchol(Y(i,i));
%         warning off MATLAB:nearlySingularMatrix
%         f = L\y(i,1); % triangular solve: L\y
%         x(i,1) = full(L'\f); % mean with original ordering
%         warning on MATLAB:nearlySingularMatrix
%         if nargout > 1
%             P = spinv(Y);
%         end

%         error('faaaaaaaaaaaaag');
%         disp('using back-slash operator');
%         warning off MATLAB:nearlySingularMatrix
%         x(i,1) = Y(i,i)\y(i,1);
%         if any(isnan(x)|~isfinite(x))
%             error('rrrrrrrr');
%         end
%         if nargout > 1
%             P(i,i) = inv(Y(i,i)); % full is used, because sparsity is already lost
%             %P = Y\speye(size(Y));
%         end
%         warning on MATLAB:nearlySingularMatrix

%    end
%end

% catch % try using spqr
%     if length(Y) > 3000;
%         error('too slow');
%     end
%     try x = spqr_solve(Y,y);
%         disp('using spqr_solve');
%         if nargout > 1
%             %P = spqr_solve(Y,I);
%             [Q,R,S] = spqr(Y);
%             %P = spqr_solve(R,Q');
%             P = R\Q';
%         end
%
%         error('too slow');
%     catch % try using ldl
%         try x = ldlsolve(Y,y);
%             disp('using ldlsolve');
%             if nargout > 1
%                 P = ldlsolve(Y,I);
%             end
%
%             %error('too slow');
%         catch % try using LU-Factorisation
%             try [L,U,p] = lu(Y,'vector');
%                 disp('using LU-Factorisation');
%                 %x = U\(L\(y(p)));
%                 f = L\y(p); % triangular solve: L\y
%                 x = full(U\f); % mean with original ordering
%                 if nargout > 1
%                     F = L\speye(size(Y));
%                     P(:,p) = (U\F);
%                 end
%
%                 error('too slow');
%             catch % try using MATLAB-PINV (SVD)
%                 try P = pinv(full(Y));
%                     disp('using pinv: PAUSED');
%                     %x = P*y;
%                     x(i) = Y(i,i)\y(i);
%
%                 catch % try using MATLAB-INV/back-slash operator
%                     disp('using back-slash operator: PAUSED');
%                     x(i) = Y(i,i)\y(i);
%                     if nargout > 1
%                         %P = inv(Y);
%                         P(i,i) = Y(i,i)\speye(size(Y));
%                     end
%
%
%                 end
%             end
%         end
%     end
% end
% end
%
%
% % mex function (using Eigen)
% i = amd(Y);
% [x, P] = mex_recover_moments(y(i), full(Y(i,i)));
% P(i,i) = P;
% x(i) = x;
%
% % sparseinv
% % Computes a subset of inv(A) for a real sparse matrix A,
% % using Takahashi's equations.
% Z, Zpattern] = sparseinv(Y);
% P = Zpattern.*Z;
% x = full(P*y);
% P(isnan(P)) = 0;
% P(isinf(P)) = 1;
