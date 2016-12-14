


% used to test if the epipoles are acceptable
% Tariq Abuhashim - August 2014, iCub

function status = check_epipoles(e1,e2,siz)


% Note: status is set to 1 when epipoles are not acceptable

% status
status = false;

% get image size
rows = siz(1);
cols = siz(2);

% % check if epipoles are in the image frames
% [FU, ~,FV]=svd(F);
% e1=FV(:,3); % left epipole e1=null(F_mapsac);
% e2=FU(:,3); % right epipole e2=null(F_mapsac');

% Put epipoles in front of camera
if e1 < 0; e1 = -e1; end
if e2 < 0; e2 = -e2; end

% non zero last coordinate
if abs(e1(3))<1e-6 && abs(e2(3))< 1e-6;
   % warning('abs(e1(3))<1e-6 && abs(e2(3))< 1e-6, epipole could be at infinity?');
    status = true;
end;

% normalise
e1 = e1/e1(3);
e2 = e2/e2(3);

% check if inside the image or non-existing
if ( e1(1) <= cols && e1(1) >= 1 && e1(2) <= rows && e1(2) >= 1 ) || ...
        ( e2(1) <= cols && e2(1) >= 1 && e2(2) <= rows && e2(2) >= 1 ) || ...
        isempty(e1) || status;
   % warning('At least one epipole is inside the image or is empty.');
    status = true;
else
    status = false;
end

