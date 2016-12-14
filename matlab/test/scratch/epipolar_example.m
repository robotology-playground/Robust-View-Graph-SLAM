%
% The Epipolar line search, similar to algorithm implemented in
% epipolargeometry.m
%

clc;
clear all;
close all;


f = 1000;
s = 0;
a = 1;
x0 = 250;
y0 = 250;

K = [f s x0; 0 a*f y0; 0 0 1];

t = [.2; .3; -.2];
phi = [2; 2.6; 2.5]*pi/180;
R = euler2DCM(phi);

d = -100:.1:100; % depth ambiguity

figure(1);
axis([0 500 0 500])
hold on;
grid on;
axis equal;

figure(2);
axis([0 500 0 500])
hold on;
grid on;
axis equal;

for itr = 1:50;
    
    epipolarline = [];
    
    figure(1);
    axis([0 500 0 500])
    [x,y] = ginput(1);
    plot(x, y, 'O');
    
    xl = [x; y; 1];
    
    for i=1:length(d)
        line = (K*R'/K)*xl + d(i)*K*R'*t;
        epipolarline(1,i) = line(1)/line(3);
        epipolarline(2,i) = line(2)/line(3);
        
        line = (K*R'/K)*xl;
        epipolar(1,i) = line(1)/line(3);
        epipolar(2,i) = line(2)/line(3);
    end
    
    figure(2);
    axis([0 500 0 500])
    plot(epipolarline(1,:), epipolarline(2,:))
    plot(epipolar(1,:), epipolar(2,:),'rO')

end


%axis([0 500 0 500])