function [imgl, imgr, eyes, neck, waist, floatingbase] = sync_images_with_encoders(cam1, cam2, hd, tor, floatingbase, freq, id)

rate = 1/freq;

disp([' - Frames captured at ',num2str(1/rate),' Hz']);
disp([' - Frames later than ',num2str(rate/2),' sec were dropped']);

if isempty(cam1)||isempty(cam2); 
    error(' - imgl or imgr data is empty'); 
end;

% right images
m1 = repmat(cam2.t,1,size(cam1.t,1));
m2 = repmat(cam1.t',size(cam2.t,1),1);
[a3,i3] = min(abs(m1-m2));
use = a3<=rate/2;
on  = find(use);
id = id(id<=length(on)); %added 26/04/2016, limits id to max images number
if ~isempty(id); on = on(id); end
imgl.image  = cam1.name(on)';
imgr.image  = cam2.name(i3(on))';
imgr.time  = a3(on);

% head encoders
if ~isempty(hd)
    m1 = repmat(hd.t,1,size(cam1.t,1));
    m2 = repmat(cam1.t',size(hd.t,1),1);
    [a1,i1] = min(abs(m1-m2));
    eyes = hd.a_eye(i1(on),:)'*pi/180;
    neck = hd.a_hd(i1(on),:)'*pi/180;
else
    eyes = [];  neck = [];
end

% waist encoders
if ~isempty(tor)
    m1 = repmat(tor.t,1,size(cam1.t,1));
    m2 = repmat(cam1.t',size(tor.t,1),1);
    [a2,i2] = min(abs(m1-m2));
    waist = tor.a(i2(on),:)'*pi/180;
else
    waist = [];
end

% floating-base
if ~isempty(floatingbase)
    m1 = repmat(floatingbase.t,1,size(cam1.t,1));
    m2 = repmat(cam1.t',size(floatingbase.t,1),1);
    [a2,i2] = min(abs(m1-m2));
    base(1:3,:) = floatingbase.dx(i2(on),:)';
    for i=1:length(on)
        R = reshape(floatingbase.R(i2(on(i)),:),3,3);
        base(4:6,i) = R2w(R)';
    end
    floatingbase = base;
else
    floatingbase = [];
end

disp([' - Number of pairs dropped = ',num2str(sum(~use))]);
disp([' - Number of pairs taken = ',num2str(length(on))]);

% old working code
%%%%%%%%%%%%%%%%%%
% function [eyes, neck, waist, imgl, imgr] = sync_images_with_encoders(cam1, cam2, hd, tor, freq, id)
% 
% rate = 1/freq;
% 
% disp([' ----> Frames captured at ',num2str(1/rate),' Hz']);
% disp([' ----> Frames later than ',num2str(rate/2),' sec were dropped']);
% 
% % head encoders
% m1 = repmat(hd.t,1,size(cam1.t,1));
% m2 = repmat(cam1.t',size(hd.t,1),1);
% [a1,i1] = min(abs(m1-m2));
% 
% % waist encoders
% m1 = repmat(tor.t,1,size(cam1.t,1));
% m2 = repmat(cam1.t',size(tor.t,1),1);
% [a2,i2] = min(abs(m1-m2));
% 
% % right images
% m1 = repmat(cam2.t,1,size(cam1.t,1));
% m2 = repmat(cam1.t',size(cam2.t,1),1);
% [a3,i3] = min(abs(m1-m2));
% 
% % drop frames with large timestamp difference
% use = a3<=rate/2;
% on  = find(use);
% if ~isempty(id); on = on(id); end
% eyes        = hd.a_eye(i1(on),:)'*pi/180;
% neck        = hd.a_hd(i1(on),:)'*pi/180;
% waist       = tor.a(i2(on),:)'*pi/180;
% imgr.image  = cam2.name(i3(on))';
% imgl.image  = cam1.name(on)';
% 
% % % take frames with specific time-stamps
% % if ~isempty(id) && length(id)<sum(on)
% %     eyes = eyes(:,id);
% %     neck = neck(:,id);
% %     waist = waist(:,id);
% %     imgl.image = imgl.image(id);
% %     imgr.image = imgr.image(id);
% % end
% 
% disp([' ----> Frames dropped ',num2str(sum(~use))]);
% disp([' ----> Frames taken ',num2str(length(on))]);
