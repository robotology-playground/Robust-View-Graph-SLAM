function [options, floatingbase] = set_images(options)
% [options, encoders] = set_images(options)
% setup folders, parameters and read images for a sequence of images.
% this is used in the patch and sequential implementations
%
% Tariq Abuhashim
% started: August 2016
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
[img,floatingbase]=load_data(path_to_img);

% sync time stamps with the left images
[cam,floatingbase]=sync_images_with_encoders(img,floatingbase,freq,image_id);

% left and right image folders from main data folder
img_folder{1} = [path_to_img,'/image/'];
img_folder{2} = [path_to_img,'/image/'];

% outputs - Do we split here? NO?
cam_left = cam;
cam_right = cam; 
options.cam_left = cam_left;
options.cam_right = cam_right;
options.img_folder = img_folder;

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
    save(strcat(options.save,'/floatingbase'),'floatingbase');
    
end

% plots
if options.verbose>1
    
    if ~isempty(floatingbase)
        figure; plot3(floatingbase(1,:),floatingbase(2,:),floatingbase(3,:)); 
        axis equal; grid on;
        figure; 
        subplot(3,1,1); plot(floatingbase(4,:)*180/pi); grid on;
        subplot(3,1,2); plot(floatingbase(5,:)*180/pi); grid on;
        subplot(3,1,3); plot(floatingbase(6,:)*180/pi); grid on;
    end
    
end
