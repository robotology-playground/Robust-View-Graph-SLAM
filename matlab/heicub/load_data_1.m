function [img,floatingbase]=load_data(folder)
%[img,floatingbase] = load_data(folder)
%
% this imports the r1 logged data, where
% 
% folder: main folder containing the data logs
% img: image data
% floatingbase: odometry data
% Tariq Abuhashim - 2016
%
% iCub

time_tx_min = [];
time_rx_min = [];
img = [];
floatingbase = [];

% read camera
if isdir([folder, '/image'])
    data = load([folder, '/image/data.log']);
    % split the images first? NO?
    if ~isempty(data);
        img.pkg = data(:, 1);
        img.tx = data(:, 2);
        img.rx = data(:, 2); % rx = tx
        names = zero_padding(data(:, 3), 8);
        img.name = strcat(names, repmat('.ppm', length(names), 1));
        img.folder = [folder, '/image/'];
        img.name = strcat(names, repmat('.ppm', length(names), 1));
        time_tx_min = [time_tx_min, img.tx(1)];
        time_rx_min = [time_rx_min, img.rx(1)];
    end
end
if isempty(img); disp('camera data is empty'); end;

% read floating base
if isdir([folder, '/odometry'])
    data = importdata([folder, '/odometry/data.log']);
    if ~isempty(data);
        floatingbase.pkg = data(:, 1);
        floatingbase.tx = data(:, 2);
        floatingbase.rx = data(:, 2); % rx = tx
        floatingbase.dx = data(:, 3:4); %translation (x,y)
        floatingbase.dx(:,3) =  0;
        floatingbase.R = floatingbase.dx*0;
        floatingbase.R(:,3) = data(:, 5)*pi/180; %rotation (heading)
        time_tx_min = [time_tx_min, floatingbase.tx(1)];
        time_rx_min = [time_rx_min, floatingbase.rx(1)];
    end
end
if isempty(floatingbase); disp('odometry data is empty'); end;

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% sync time stamps to first left image and remove negative times
if 1;
    % use tx time
    time_tx_min = 0; % uncomment this line to set reference time to start at 0
    if ~isempty(img);
        img.t = img.tx - time_tx_min;
        img.t = img.t(img.t>0,:);
        img.pkg = img.pkg(img.t>0,:);
        img.tx = img.tx(img.t>0,:);
        img.rx = img.rx(img.t>0,:);
        img.name = img.name(img.t>0,:);
    end
    if ~isempty(floatingbase);
        floatingbase.t = floatingbase.tx - time_tx_min;
        floatingbase.t = floatingbase.t(floatingbase.t>0,:);
        floatingbase.pkg = floatingbase.pkg(floatingbase.t>0,:);
        floatingbase.tx = floatingbase.tx(floatingbase.t>0,:);
        floatingbase.rx = floatingbase.rx(floatingbase.t>0,:);
        floatingbase.dx = floatingbase.dx(floatingbase.t>0,:);
        floatingbase.R = floatingbase.R(floatingbase.t>0,:);
    end
else
    % use rx time
    time_rx_min = 0; % comment this line to set reference time to start at 0
    if ~isempty(img);
        img.t = img.rx - time_rx_min;
        img.t = img.t(imgl.t>0,:);
        img.pkg = img.pkg(imgl.t>0,:);
        img.tx = img.tx(imgl.t>0,:);
        img.rx = img.rx(imgl.t>0,:);
        img.name = img.name(imgl.t>0,:);
    end
    if ~isempty(floatingbase);
        floatingbase.t = floatingbase.rx - time_rx_min;
        floatingbase.t = floatingbase.t(floatingbase.t>0,:);
        floatingbase.pkg = floatingbase.pkg(floatingbase.t>0,:);
        floatingbase.tx = floatingbase.tx(floatingbase.t>0,:);
        floatingbase.rx = floatingbase.rx(floatingbase.t>0,:);
        floatingbase.dx = floatingbase.dx(floatingbase.t>0,:);
        floatingbase.R = floatingbase.R(floatingbase.t>0,:);
    end
end

function str = zero_padding(number, strlength)

% pads the array of numbers "number" by zeros such that all elements have
% the same length.
%
% Tariq Abuhashim, 2013

if size(number,2) > size(number, 1);
    number = number';
end
if size(number, 2) ~= 1;
    error('number must be a column vector');
end
if any( ~isreal(number) | number < 0 | round(number) ~= number );
    error('number must be positive integer');
end

str = cell(size(number, 1), 1);
for i = 1:size(number, 1);
    %str{i} = sprintf('%d', number(i));
    str{i} = num2str(number(i));
    difflength = strlength - length(str{i});

    for j = 1:difflength
        %str{i} = strcat('0', str{i});
        str{i} = ['0', str{i}];
    end
end
