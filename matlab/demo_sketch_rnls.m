clear;

% a simulation for robust nonlinear least-squares
%y = @(x) 2*(exp(x.*x - 0.4) - 1.2*x.*x - .2*x); % function
%H = @(x) 2*(2.*x.*exp(x.*x - 0.4) - 2.4*x - .2); % jacobian
y = @(x) - x.^5 + 8*x.^3 - 15*x + 0; % function
H = @(x) - 5.*x.^4 + 24*x.^2 - 15 + 0; % jacobian
close all

xf = -.5:.01:2.2; % function
yf = feval(y,xf); % function

figure;
plot(xf,yf,'b','linewidth',2); % function
axis([-.5 2.2 -20 10])
set(gca,'XTickLabel',[],'YTickLabel',[]);
hold on;
box on;

% xhat
xhat = 0.70; % estimate
xn = linspace(xhat-.3,xhat+.3,length(xf)); % distribution - x
xnorm = normpdf(xn,xhat,.1);
plot(xn,xnorm-20,'b','linewidth',2);

% xs
xs = xhat;
ys = feval(y,xs);
J = feval(H,xs);
yhat = J*(xhat-xs) + ys;
plot(xhat,yhat,'bO','linewidth',2);
plot([xs xs],[-20 ys],'k--');

% tangent
xl = linspace(xs-.4,xs+.4,10);
yl = J*(xl-xs) + ys;
plot(xl,yl,'r','linewidth',2);

% measurement
yn_o = J*(xn-xhat) + yhat; % distribution - y
ynorm_o = abs(J*xnorm); % distribution - y
plot(ynorm_o/20-.5,yn_o,'r','linewidth',2); % distribution - y
plot([xhat xhat],[-20 yhat],'k--');
plot([-.5 xhat],[yhat yhat],'k--');

%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
plot(xf,yf,'b','linewidth',2); % function
axis([-.5 2.2 -20 10])
set(gca,'XTickLabel',[],'YTickLabel',[]);
hold on;
box on;

% xhat
plot(xn,xnorm-20,'b','linewidth',2);

% linearisation point
xs = 1.15;
ys = feval(y,xs);
J = feval(H,xs);
plot(xs,ys,'bO','linewidth',2);
plot([xs xs],[-20 ys],'k--');

% previous measurement
plot(ynorm_o/20-.5,yn_o,'r--','linewidth',2); % distribution - y

% tangent
xl = linspace(xs-.7,xs+.3,10);
yl = J*(xl-xs) + ys;
plot(xl,yl,'r','linewidth',2);

% measurement
yhat = J*(xhat-xs) + ys;
plot(xhat,yhat,'bO','linewidth',2);
yn_o = J*(xn-xhat) + yhat; % distribution - y
ynorm_o = abs(J*xnorm); % distribution - y
plot(ynorm_o/50-.5,yn_o,'r','linewidth',2); % distribution - y
plot([xhat xhat],[-20 yhat],'k--');
plot([-.5 xhat],[yhat yhat],'k--');


return % old code

xs_1 = xhat; % red linearisation point (Kalman filter)
ys_1 = feval(y,xs_1);
J_1 = feval(H,xs_1);

plot(xs_1,ys_1,'rO','linewidth',2);
plot([xs_1 xs_1],[0 ys_1],'k--');

xs_2 = 0.4; % green linearisation point (robust nonlinear least squares)
ys_2 =  feval(y,xs_2);
J_2 = feval(H,xs_2);

plot(xs_2,ys_2,'gO','linewidth',2);
plot([xs_2 xs_2],[0 ys_2],'k--');

xl = linspace(0,1,10); % tangent lines
yl_1 = J_1*(xl-xs_1) + ys_1;
yl_2 = J_2*(xl-xs_2) + ys_2;

plot(xl,yl_1,'r','linewidth',2)
plot(xl,yl_2,'g','linewidth',2)

yn_1 = J_1*(xn-xs_1) + ys_1; % distribution - y
ynorm_1 = abs(J_1*xnorm);
yn_2 = J_2*(xn-xs_2) + ys_2;
ynorm_2 = abs(J_2*xnorm);

plot(ynorm_1/50,yn_1,'r','linewidth',2); % red distribution
plot(ynorm_2/50,yn_2,'g','linewidth',2); % green distribution

yhat_1 = J_1*(xhat-xs_1) + ys_1;
yhat_2 = J_2*(xhat-xs_2) + ys_2;

plot([0 xhat],[yhat_1 yhat_1],'k--');
plot([0 xhat],[yhat_2 yhat_2],'k--');