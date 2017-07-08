function [options, encoders, floatingbase] = set_images(options)
% [options, encoders] = set_images(options)
% setup folders, parameters and read images for a sequence of images.
% this is used in the patch and sequential implementations
%
% Tariq Abuhashim
% started: September 2014
%
% iCub - Koroibot

if isstruct(options)
    if isfield(options,'folder') && isfield(options,'freq')
        path_to_img = options.folder;
        freq = options.freq;
        if isfield(options,'first_image');   
            start = options.first_image;
        else
            start = 0;                        
        end;
        if isfield(options,'last_image');    
            stop = options.last_image;
        else
            stop = 0;                        
        end;
        if isfield(options,'steps');         
            steps = options.steps;
        else
            steps = 0;                        
        end;
        image_id = start:steps:stop; % images used
    end
else
    error('Input should be either a structure with at least two fields (folder, freq)');
end

% get icub data
[imu,hd,lh,rh,la,ra,tor,imgl,imgr,floatingbase]=load_data(path_to_img);

% sync time stamps with the left images
[cam_left,cam_right,eyes,neck,waist,floatingbase]=sync_images_with_encoders(imgl,imgr,hd,tor,floatingbase,freq,image_id);

% left and right image folders from main data folder
img_folder{1} = [path_to_img,'/img/left/'];
img_folder{2} = [path_to_img,'/img/right/'];

% outputs
options.cam_left = cam_left;
options.cam_right = cam_right;
options.img_folder = img_folder;
encoders = [eyes; neck; waist];

% save to .mat file:
if isfield(options,'save');
    % ============= % create save directory
    if isfield(options,'save');
        if ~isdir(strcat(options.save));
            mkdir(strcat(options.save));
        end
    end
    save(strcat(options.save,'/options'),'options');
    save(strcat(options.save,'/cams'),'cam_left','cam_right');
    save(strcat(options.save,'/encoders'),'eyes', 'neck', 'waist');
    save(strcat(options.save,'/floatingbase'),'floatingbase');
end

% plots
if options.verbose>1
    if ~isempty(eyes)&&~isempty(neck)
        figure;
        plot([eyes; neck]'*180/pi,'linewidth', 2);grid on;
        %xlabel('images','interpreter','latex','fontsize',20)
        %ylabel('angles (deg.)','interpreter','latex','fontsize',20)
        %set(gca,'fontsize',16);
        legend('tilt','pan','vergence','head roll','head pitch','head yaw');
    end
    if ~isempty(eyes)&&~isempty(neck)
        figure;
        plot(waist'*180/pi,'linewidth', 2);grid on;
        legend('waist roll','waist pitch','waist yaw');
        %xlabel('images'); ylabel('angles in degrees');
    end
    if ~isempty(floatingbase)
        figure; plot3(floatingbase(1,:),floatingbase(2,:),floatingbase(3,:)); 
        axis equal; grid on;
        figure; 
        subplot(3,1,1); plot(floatingbase(4,:)*180/pi); grid on;
        subplot(3,1,2); plot(floatingbase(5,:)*180/pi); grid on;
        subplot(3,1,3); plot(floatingbase(6,:)*180/pi); grid on;
    end
    figure;
    t = min([imgr.t;imgl.t]);
    plot(imgl.t-t);hold on;plot(imgr.t-t,'r');
    figure;
    subplot(2,1,1);plot(diff(imgl.t));axis tight;
    subplot(2,1,2);plot(diff(imgr.t),'r');axis tight;
    figure;
    plot(diff(hd.t));hold on;plot(diff(tor.t),'r');
    title('encoders time step');
    drawnow;
    figure;
    plot(cam_right.time);hold on;
    title('Right to left time step');
    drawnow;
end
