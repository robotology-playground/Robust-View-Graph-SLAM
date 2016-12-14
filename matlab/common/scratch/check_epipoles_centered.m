


% used to test if the epipoles are acceptable for the case of centred image
% coordinates
% Tariq Abuhashim - August 2014, iCub

function status=check_epipoles_centered(e1,e2,siz)


% get image size
rows=siz(1);
cols=siz(2);

% % check if epipoles are in the image frames
% [FU, ~,FV]=svd(F);
% e1=FV(:,3); % left epipole e1=null(F_mapsac);
% e2=FU(:,3); % right epipole e2=null(F_mapsac');

% Put epipoles in front of camera
if e1 < 0; e1 = -e1; end
if e2 < 0; e2 = -e2; end

% non zero last coordinate
if abs(e1(3))<1e-6 && abs(e2(3))< 1e-6;
    warning('abs(e1(3))<1e-6 && abs(e2(3))< 1e-6, epipole could be at infinity ? \n')
end

% normalise
e1 = e1/e1(3);
e2 = e2/e2(3);

% check if inside the image or non-existing
nx = cols/2;  ny = rows/2;
if ( (e1(1)>=-nx && e1(1)<=nx && e1(2)>=-ny && e1(2)<=ny) || ...
        ( e2(1)>=-nx && e2(1)<=nx && e2(2)>=-ny && e2(2)<=ny) ) || ...
        isempty(e1);
    status=1;
else
    status=0;
end

