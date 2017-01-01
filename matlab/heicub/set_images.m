function [options,encoders,floatingbase]=set_images(options)
%[options,encoders,floatingbase]=set_images(options)
%
% setup folders, parameters and read images for a sequence of images.
% this is used in the patch and sequential implementations
%
%	For iCub@heidelberg.
%
% Tariq Abuhashim
% started: September 2014
% t.abuhashim@gmail.com
%
% iCub - Koroibot

if isstruct(options)
    if isfield(options,'folder') && isfield(options,'freq')
        path_to_img = options.folder;
        freq = options.freq;
	start = 0;
        if isfield(options,'first_image')
            start = options.first_image;                      
        end
	stop = 0;
        if isfield(options,'last_image')
            stop = options.last_image;                       
        end
	steps = 0; 
        if isfield(options,'steps')
            steps = options.steps;                     
        end
        image_id = start:steps:stop; % images used
    end
else
    error('Input should be a structure with at least two fields (folder, freq)');
end

% get icub data
[tor,imgl,imgr,floatingbase]=load_data(path_to_img);
    
% sync time stamps with the left images
[cam_left,cam_right,waist,floatingbase]=sync_images_with_encoders(imgl,imgr,tor,floatingbase,freq,image_id);

% left and right image folders from main data folder
img_folder{1} = [path_to_img,'/img/left/'];
img_folder{2} = [path_to_img,'/img/right/'];

% outputs
options.cam_left = cam_left;
options.cam_right = cam_right;
options.img_folder = img_folder;
encoders([7 8 9],:) = waist;

% save to .mat file:
if isfield(options,'save')
    % ============= % create save directory
    if isfield(options,'save')
        if ~isdir(strcat(options.save))
            mkdir(strcat(options.save));
        end
    end
    save(strcat(options.save,'/options'),'options');
    save(strcat(options.save,'/cams'),'cam_left','cam_right');
    save(strcat(options.save,'/encoders'), 'waist');
    save(strcat(options.save,'/floatingbase'),'floatingbase');
end

% plots
if options.verbose>1
    
    if ~isempty(waist)
        figure;
        plot(waist'*180/pi,'linewidth', 2);grid on;
        legend('waist roll','waist pitch','waist yaw');
        %xlabel('images'); ylabel('angles in degrees');
    end
    
    if ~isempty(floatingbase)
        figure; plot3(floatingbase(1,:),floatingbase(2,:),floatingbase(3,:)); 
        axis equal; grid on;
        title('floating position')
        figure;
        subplot(3,1,1); plot(floatingbase(4,:)); grid on;
        title('floating orientation')
        subplot(3,1,2); plot(floatingbase(5,:)); grid on;
        subplot(3,1,3); plot(floatingbase(6,:)); grid on;
    end
    
    figure;
    plot(cam_left.time);hold on;plot(cam_right.time,'r');
    title('left and right camera time stamps');
    
    figure;
    subplot(2,1,1);plot(diff(imgl.time));axis tight;
    title('left camera time step');
    subplot(2,1,2);plot(diff(imgr.time),'r');axis tight;
    title('right camera time step');
    
    figure;
    plot(diff(tor.time),'r');
    title('encoders time step');
    drawnow;
    
    figure;
    plot(cam_right.time);hold on;
    title('Right to left time step');
    drawnow;
    
end
