function [imgl, imgr] = load_data(folder)
%[imgl, imgr] = load_data(folder)
%
% Imports the HRP2 logged data
%
%	For HRP2@Toulouse.
%
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Auckland, New Zealand, 2017

imgl = [];
imgr = [];

% read left camera
if isdir([folder, '/left'])
    filenames = dir([folder,'/left/*.pgm']);
    if ~isempty(filenames)
    	imgl.name = strings(1,size(filenames,1));
    	for i=1:size(filenames,1); 
        	imgl.name(i) = filenames(i).name;
        end
        imgl.folder = [folder,'/left/'];
    end
end
if isempty(imgl)
    disp('left camera data is empty');
end

% read right camera
if isdir([folder, '/right'])
    filenames = dir([folder,'/right/*.pgm']);
    if ~isempty(filenames)
        imgr.name = strings(1,size(filenames,1));
    	for i=1:size(filenames,1); 
        	imgr.name(i) = filenames(i).name;
        end
        imgr.folder = [folder,'/right/'];
    end
end
if isempty(imgr)
    disp('right camera data is empty');
end
