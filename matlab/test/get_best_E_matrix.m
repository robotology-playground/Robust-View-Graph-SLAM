
function [E, J, C, tr] = get_best_E_matrix(p1, p2, sigma_r)

sigma_r = sigma_r/410;

J = [];
for i = 1:10;
    Ji = numerical_jacobian_i(@test_Evec,[],i,[],...
        p1(1:2,1),p1(1:2,2),p1(1:2,3),p1(1:2,4),p1(1:2,5),...
        p2(1:2,1),p2(1:2,2),p2(1:2,3),p2(1:2,4),p2(1:2,5));
    J = [J Ji];
end

J

R = sigma_r*eye(20);
tr = zeros(1,size(J,1)/9);
for i = 1:size(J,1)/9
    idx = (i-1)*9 + (1:9);
    C = J(idx,:)*R*J(idx,:)';
    tr(i) = trace(C);
end

plot(tr); pause

[~, i] = min(tr);
idx = (i-1)*9 + (1:9);
Evec = calibrated_fivepoint(p1, p2);
E = reshape(Evec(:,i), 3, 3);
J = J(idx,:);
C = J*R*J';
tr = tr(i);

function out = test_Evec(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9)
p1_sample = [p0 p1 p2 p3 p4]; p1_sample = pextend(p1_sample);
p2_sample = [p5 p6 p7 p8 p9]; p2_sample = pextend(p2_sample);
Evec = calibrated_fivepoint(p1_sample, p2_sample);
out = Evec(:);