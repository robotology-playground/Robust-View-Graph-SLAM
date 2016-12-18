function [tor, imgl, imgr, floatingbase] = load_data(folder)
%[tor, imgl, imgr] = load_data(folder)
%
% Imports the HeiCub logged data
%
% iCub Facility
% Tariq Abuhashim - 2014
%

time_min = +Inf;
tor  = [];
imgl = [];
imgr = [];
floatingbase = [];

% read left camera
if isdir([folder, '/img/left'])
    data = load([folder, '/img/left/data.log']);
    if ~isempty(data)
        imgl.pkg = data(:,1);
        imgl.time = data(:,2);
        names = zero_padding(data(:,3),8);
        imgl.name = strcat(names,repmat('.ppm',length(names),1));
        imgl.folder = [folder,'/img/left/'];
    end
end
if isempty(imgl)
    disp('left camera data is empty');
else % remove unreadable files
    idx = [];
    count = 0;
    for i = 1:length(imgl.name)
        if exist(strcat([folder, '/img/left'],'/',imgl.name{i}),'file')
            idx = [idx; i];
        else
            count = count+1;
        end
    end
    if count>0
        fprintf('%d left images do not exist\n',count);
    end
    imgl.pkg = imgl.pkg(idx);
    imgl.time = imgl.time(idx);
    imgl.name = imgl.name(idx);
    time_min = min(time_min,imgl.time(1));
end

% read right camera
if isdir([folder, '/img/right'])
    data = load([folder, '/img/right/data.log']);
    if ~isempty(data)
        imgr.pkg = data(:, 1);
        imgr.time = data(:, 2);
        names = zero_padding(data(:, 3), 8);
        imgr.name = strcat(names, repmat('.ppm', length(names), 1));
        imgr.folder = [folder, '/img/right/'];
    end
end
if isempty(imgr)
    disp('right camera data is empty');
else % remove unreadable files
    idx = [];
    count = 0;
    for i = 1:length(imgr.name)
        if exist(strcat([folder, '/img/right'],'/',imgr.name{i}),'file')
            idx = [idx; i];
        else
            count = count+1;
        end
    end
    if count>0
        fprintf('%d right images do not exist\n',count);
    end
    imgr.pkg = imgr.pkg(idx);
    imgr.time = imgr.time(idx);
    imgr.name = imgr.name(idx);
    time_min = min(time_min,imgr.time(1));
end

% read torso
if isdir([folder, '/torso'])
    data = importdata([folder, '/torso/data.log']);
    if ~isempty(data)
        tor.pkg = data(:, 1);
        tor.time = data(:, 2);
        tor.a = data(:, [4 5 3]); % body orientation (data = [4-yaw 5-roll 6-pitch], a = [roll pitch yew];
        time_min = min(time_min,tor.time(1));
    end
end
if isempty(tor); disp('torso data is empty'); end;

% read floating base
if isdir([folder, '/floatingBase'])
    data = importdata([folder, '/floatingBase/data.log']);
    if ~isempty(data)
        floatingbase.pkg = data(:, 1);
        floatingbase.time = data(:, 2);
	floatingbase.translation = data(:, 6:8)/1000; %translation vector
        floatingbase.rotation = data(:, 3:5)*pi/180; %rotation vector
        time_min = min(time_min,floatingbase.time(1));
    end
end
if isempty(floatingbase); disp('floatingbase data is empty'); end;

% sync time stamps to first left image and remove negative times
% use tx time
% time_min = 0; % uncomment this line to set reference time to start at 0
if ~isempty(imgl)
    t = imgl.time - time_min;
    imgl.time = imgl.time(t>0,:);
    imgl.pkg = imgl.pkg(t>0,:);
    imgl.name = imgl.name(t>0,:);
end
if ~isempty(imgr)
    t = imgr.time - time_min;
    imgr.time = imgr.time(t>0,:);
    imgr.pkg = imgr.pkg(t>0,:);
    imgr.name = imgr.name(t>0,:);
end
if ~isempty(tor)
    t = tor.time - time_min;
    tor.time = tor.time(t>0,:);
    tor.pkg = tor.pkg(t>0,:);
    tor.a = tor.a(t>0,:);
end
if ~isempty(floatingbase)
    t = floatingbase.time - time_min;
    floatingbase.pkg = floatingbase.pkg(t>0,:);
    floatingbase.time = floatingbase.time(t>0,:);
    floatingbase.translation = floatingbase.translation(t>0,:);
    floatingbase.rotation = floatingbase.rotation(t>0,:);
end

function str = zero_padding(number, strlength)

% Pad the array of numbers "number" by zeros such that all elements have
% the same length.
%
% Tariq Abuhashim, 2013

if size(number,2) > size(number, 1)
    number = number';
end
if size(number, 2) ~= 1
    error('number must be a column vector');
end
if any( ~isreal(number) | number < 0 | round(number) ~= number )
    error('number must be positive integer');
end

str = cell(size(number, 1), 1);
for i = 1:size(number, 1)
    %str{i} = sprintf('%d', number(i));
    str{i} = num2str(number(i));
    difflength = strlength - length(str{i});
    
    for j = 1:difflength
        %str{i} = strcat('0', str{i});
        str{i} = ['0', str{i}];
    end
end
