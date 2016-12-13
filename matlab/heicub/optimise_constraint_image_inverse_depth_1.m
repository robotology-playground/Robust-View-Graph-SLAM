function [y,Y,sw,x]=optimise_constraint_image_inverse_depth(Ct,C,sw,xs,options,ncams)
%[y,Y,sw,x]=optimise_constraint_image_inverse_depth(Ct,C,sw,xs,options)
%
% Tariq Abuhashim, 2016.
% iCub - Koroibot

% verbose
if 1;%options.verbose > 0
    %disp('  ')
    disp(['Trusted constraints: ' num2str(length(Ct))])
    disp(['ON constraints: ' num2str(sum(sw == 1))])
    disp(['OFF constraints: ' num2str(sum(sw == 0))])
end

% Initialise estimate with switched ON constraints
[y, Y] = initialise_info_matrix(Ct, xs, ncams); % sw is for C, not for Ct

% Generate image constraints information
% H = sparse(2*(length(xs)-6*ncams),length(xs));
% for k = 1:length(C)
%     i = 6*ncams+C(k).kpt;
%     c = (C(k).cam-1)*6;
%     h = mex_observation_model_jacobian_inverse_depth_Mviews( xs, C(k).p1, i, c+1 );
%     H(2*C(k).kpt-1,[c+(1:6) i]) = h(1,:);
%     H(2*C(k).kpt,[c+(1:6) i]) = h(1,:);
% end
% clf; spy(H); pause;
if 0
    C = generate_image_constraints_info_inverse_depth(C, xs, ncams);
else
    C = mex_generate_constraints_info_Mviews(C, xs, ncams);
end

% Include measurements that are initially ON, but not trusted
npts = length(xs) - 6*ncams;
if 0
    [yon, Yon] = update_info_matrix_inverse_depth(C(sw == 1), npts, ncams);
else
    [yon, Yon] = mex_update_info_matrix_Mviews(C(sw == 1), npts, ncams);
end
y = y + yon;
Y = Y + Yon;

% print to file
%x = recover_moments(y, Y);
%for i = 1:6*ncams; fprintf(options.fid, '%.4f ', x(i)); end
%fprintf(options.fid, [num2str(length(sw)) ' ' num2str(sum(sw)) '\n']);

% verbose (show map uncertainty)
%if options.verbose > 1;
%    plot_scan_with_uncertainty(y, Y, options.p, ncams, 'initial');
%    pause;
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FIXME: Both mex functions for information addition and subtraction has an%
%issue with matrix inversion code (mex_recover_moments):                  %
%   1. It only supports the positive-definite case.                       %
%   2. Lines 280-283 has a memory allocation (assignment) error.          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply remaining constraints with residual switching
tic;
if 1 % Add
    [y, Y, sw] = constraints_addition_inverse_depth(y, Y, C, sw, xs, options, ncams);
else
    [y, Y, sw] = mex_constraints_addition_inverse_depth_Mviews(y, Y, C, sw, xs, ncams);
end
while any(sw==1)
    if 1 % Subtract
        [y, Y, sw, ~, converged] = constraints_removal_inverse_depth(y, Y, C, Ct, sw, xs, options, ncams);
    else
    	[y, Y, sw] = mex_constraints_subtraction_inverse_depth_Mviews(y, Y, C, sw, xs, ncams);
    end
    if converged; break; end
end
if 1 % Solve
    x = recover_moments(y, Y);
else
    x = mex_recover_moments(y, Y);
end
toc

%for i = 1 : ncams
%   fprintf('%.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n', x((i-1)*6+(1:3))', x((i-1)*6+(4:6))'*180/pi);
%end
%pause

% print to file
%for i = 1:6*ncams; fprintf(options.fid, '%.4f ', x(i)); end;
%fprintf(options.fid, [num2str(length(sw)) ' ' num2str(sum(sw)) '\n']);

% verbose (show map uncertainty)
%if options.verbose > 1;
%    plot_scan_with_uncertainty(y, Y, options.p, ncams, 'subtract');
%    pause;
%end

%
%
function plot_scan_with_uncertainty(y, Y, p, ncams, filename)
N = length(Y)-6*ncams;
[x, P] = recover_moments(y, Y);
r = 1./x(6*ncams+1:size(x,1))';
s = sqrt(diag(P((1:N)+6*ncams,(1:N)+6*ncams)))';
p0 = get_scan_from_range(p,r);
p1 = get_scan_from_range(p,r-s);
p2 = get_scan_from_range(p,r+s);
clf; plot(p0(1,r>0),p0(3,r>0),'b+','markersize',2); hold on;
plot([p1(1,r>0);p2(1,r>0)],[p1(3,r>0);p2(3,r>0)],'m');
axis equal; axis([-1 1 0 2]); grid on; box on;
set(gca, 'fontsize', 14);
xlabel('x-direction, mm', 'fontsize', 14);
ylabel('z-direction, mm', 'fontsize', 14);
set(gcf, 'Color', 'w');
export_fig([filename,'.png']);
%
%
function p = get_scan_from_range(p, r)
if size(p,1) < 3; p = pextend(p); end;
rim = sqrt(sum(p(1:2,:).^2) + 1);
d = r./rim;
p(1,:) = p(1,:).*d;
p(2,:) = p(2,:).*d;
p(3,:) = p(3,:).*d;
