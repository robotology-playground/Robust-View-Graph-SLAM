function [img, base] = sync_images_with_encoders(cam, floatingbase, freq, id)

rate = 1/freq;
disp([' ----> Frames captured at ',num2str(1/rate),' Hz']);
disp([' ----> Frames with odometry data later than ',num2str(rate/2),' sec were dropped']);

if isempty(cam); 
    error('img data is empty'); 
end
if isempty(floatingbase); 
    error('floatingbase data is empty'); 
end

m1 = repmat(floatingbase.t,1,size(cam.t,1));
m2 = repmat(cam.t',size(floatingbase.t,1),1);
[a2,i2] = min(abs(m1-m2));
use = a2<=rate/2;
on = find(use);
id = id(id<=length(on)); % added 26/04/2016, limits id to max available number of images (if id>length(on))
if ~isempty(id); on = on(id); end
img.image  = cam.name(on)';
base(1:3,:) = floatingbase.dx(i2(on),:)';
base(4:6,:) = floatingbase.R(i2(on),:)';

disp([' ----> Frames dropped ',num2str(sum(~use))]);
disp([' ----> Frames taken ',num2str(length(on))]);
