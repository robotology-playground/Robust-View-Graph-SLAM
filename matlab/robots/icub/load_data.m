function [imu, hd, lh, rh, la, ra, tor, imgl, imgr, floatingbase] = load_data(folder)
%[imu, hd, lh, rh, la, ra, tor, imgl, imgr] = load_data(folder)
%
% this importdatas the iCub logged data
%
% Tariq Abuhashim - 2014
%
% iCub

time_tx_min = [];
time_rx_min = [];
imu  = [];
hd   = [];
lh   = [];
rh   = [];
la   = [];
ra   = [];
tor  = [];
imgl = [];
imgr = [];
floatingbase = [];

% read left camera
if isdir([folder, '/img/left'])
    data = load([folder, '/img/left/data.log']);
    if ~isempty(data);
        imgl.pkg = data(:, 1);
        imgl.tx = data(:, 2);
        imgl.rx = data(:, 3);
        names = zero_padding(data(:, 4), 8);
        imgl.name = strcat(names, repmat('.ppm', length(names), 1));
        imgl.folder = [folder, '/img/left/'];
        time_tx_min = [time_tx_min, imgl.tx(1)];
        time_rx_min = [time_rx_min, imgl.rx(1)];
    end
end
if isempty(imgl); 
    disp(' - left camera data is empty'); 
else % remove unreadable files
    idx = [];
    count = 0;
    for i = 1:length(imgl.name)
        if exist(strcat([folder, '/img/left'],'/',imgl.name{i}),'file')
            idx = [idx; i];
        else
           %disp('image removed');
           count = count+1;
        end
    end
    if count> 0
    	fprintf(' - %d left images do not exist',count);
    end
    imgl.pkg = imgl.pkg(idx);
    imgl.tx = imgl.tx(idx);
    imgl.rx = imgl.rx(idx);
    imgl.name = imgl.name(idx);
    time_tx_min = [time_tx_min, imgl.tx(1)];
    time_rx_min = [time_rx_min, imgl.rx(1)];
end

% read right camera
if isdir([folder, '/img/right'])
    data = load([folder, '/img/right/data.log']);
    if ~isempty(data);
        imgr.pkg = data(:, 1);
        imgr.tx = data(:, 2);
        imgr.rx = data(:, 3);
        names = zero_padding(data(:, 4), 8);
        imgr.name = strcat(names, repmat('.ppm', length(names), 1));
        imgr.folder = [folder, '/img/right/'];
        time_tx_min = [time_tx_min, imgr.tx(1)];
        time_rx_min = [time_rx_min, imgr.rx(1)];
    end
end
if isempty(imgr); 
    disp(' - right camera data is empty'); 
else % remove unreadable files
    idx = [];
    count = 0;
    for i = 1:length(imgr.name)
        if exist(strcat([folder, '/img/right'],'/',imgr.name{i}),'file')
            idx = [idx; i];
        else
           %disp('image removed');
           count = count+1;
        end
    end
    if count> 0
    	fprintf(' - %d right images do not exist',count);
    end
    imgr.pkg = imgr.pkg(idx);
    imgr.tx = imgr.tx(idx);
    imgr.rx = imgr.rx(idx);
    imgr.name = imgr.name(idx);
    time_tx_min = [time_tx_min, imgr.tx(1)];
    time_rx_min = [time_rx_min, imgr.rx(1)];
end


% read inertial
if isdir([folder, '/inertial'])
    data = importdata([folder, '/inertial/data.log']);
    if ~isempty(data);
        imu.pkg = data(:, 1);
        imu.tx = data(:, 2);
        imu.rx = data(:, 3);
        imu.a = data(:, 4:6); % euler angles [3]: deg
        imu.f = data(:, 7:9); % linear acceleration [3]: m/s^2
        imu.w = data(:, 10:12); % angular speed [3]: deg/s
        imu.B = data(:, 13:15); % magnetic field [3]: arbitrary units
        time_tx_min = [time_tx_min, imu.tx(1)];
        time_rx_min = [time_rx_min, imu.rx(1)];
    end
end
if isempty(imu); disp(' - imu data is empty'); end;

% read head
if isdir([folder, '/head'])
    data = importdata([folder, '/head/data.log']);
    if ~isempty(data);
        hd.pkg = data(:, 1);
        hd.tx = data(:, 2);
        hd.rx = data(:, 3);
        hd.a_hd = data(:, [5 4 6]); % euler angles, data = [4-pitch 5-roll 6-yaw], a_hd = [roll pitch yew];
        hd.a_eye = data(:, 7:9); % eye orientation, data = [7-tile, 8-pan (vs), 9-vergence (vg)]
        time_tx_min = [time_tx_min, hd.tx(1)];
        time_rx_min = [time_rx_min, hd.rx(1)];
    end
end
if isempty(hd); disp(' - head data is empty'); end;

% read torso
if isdir([folder, '/torso'])
    data = importdata([folder, '/torso/data.log']);
    if ~isempty(data);
        tor.pkg = data(:, 1);
        tor.tx = data(:, 2);
        tor.rx = data(:, 3);
        tor.a = data(:, [5 6 4]); % body orientation (data = [4-yaw 5-roll 6-pitch], a = [roll pitch yew];
        time_tx_min = [time_tx_min, tor.tx(1)];
        time_rx_min = [time_rx_min, tor.rx(1)];
    end
end
if isempty(hd); disp(' - torso data is empty'); end;

% read floating base
if isdir([folder, '/floating_base'])
    data = importdata([folder, '/floating_base/data.log']);
    if ~isempty(data);
        floatingbase.pkg = data(:, 1);
        floatingbase.tx = data(:, 3);
        floatingbase.rx = data(:, 3);
        floatingbase.dx = data(:, 4:6); %translation vector
        floatingbase.R = data(:, 7:15); %rotation vector
        time_tx_min = [time_tx_min, floatingbase.tx(1)];
        time_rx_min = [time_rx_min, floatingbase.rx(1)];
    end
end
if isempty(floatingbase); disp(' - floatingbase data is empty'); end;

% % read left hand
% data = importdata([folder, '/cartesian/left/data.log']);
% if isempty(data);
%     lh = [];
%     disp('left hand data is empty');
% else
%     lh.pkg = data(:, 1);
%     lh.tx = data(:, 2);
%     lh.rx = data(:, 3);
%     lh.x = data(:, 4:6); % xyz location
%     lh.v = data(:, 7:9); % orientation vector
%     lh.mag = data(:, 10); % orientation vector magnitude
%     time_tx_min = [time_tx_min, lh.tx(1)];
%     time_rx_min = [time_rx_min, lh.rx(1)];
% end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% minimum time

%time_tx_min = min(time_tx_min);
%time_tx_min = imgl.tx(30); % removes some of the starting images
%time_tx_min = 0;

%time_tx_min = min([imu.tx(1), hd.tx(1), lh.tx(1),...
%    rh.tx(1), la.tx(1), ra.tx(1), tor.tx(1), imgl.tx(1), imgr.tx(1)]);
%time_tx_max = min([imu.tx(end), hd.tx(end), lh.tx(end),...
%    rh.tx(end), la.tx(end), ra.tx(end), tor.tx(end), imgl.tx(end), imgr.tx(end)]);
%time_rx_min = min([imu.rx(1), hd.rx(1), lh.rx(1),...
%    rh.rx(1), la.rx(1), ra.rx(1), tor.rx(1), imgl.rx(1), imgr.rx(1)]);
%time_rx_max = min([imu.rx(end), hd.rx(end), lh.rx(end),...
%    rh.rx(end), la.rx(end), ra.rx(end), tor.rx(end), imgl.rx(end), imgr.rx(end)]);

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% sync time stamps to first left image and remove negative times
if 1;
    % use tx time
    time_tx_min = 0; % uncomment this line to set reference time to start at 0
    if ~isempty(imgl);
        imgl.t = imgl.tx - time_tx_min;
        imgl.t = imgl.t(imgl.t>0,:);
        imgl.pkg = imgl.pkg(imgl.t>0,:);
        imgl.tx = imgl.tx(imgl.t>0,:);
        imgl.rx = imgl.rx(imgl.t>0,:);
        imgl.name = imgl.name(imgl.t>0,:);
    end
    if ~isempty(imgr);
        imgr.t = imgr.tx - time_tx_min;
        imgr.t = imgr.t(imgr.t>0,:);
        imgr.pkg = imgr.pkg(imgr.t>0,:);
        imgr.tx = imgr.tx(imgr.t>0,:);
        imgr.rx = imgr.rx(imgr.t>0,:);
        imgr.name = imgr.name(imgr.t>0,:);
    end
    if ~isempty(imu);
        imu.t = imu.tx - time_tx_min;
        imu.t = imu.t(imu.t>0,:);
        imu.pkg = imu.pkg(imu.t>0,:);
        imu.tx = imu.tx(imu.t>0,:);
        imu.rx = imu.rx(imu.t>0,:);
        imu.a = imu.a(imu.t>0,:);
        imu.f = imu.f(imu.t>0,:);
        imu.w = imu.w(imu.t>0,:);
        imu.B = imu.B(imu.t>0,:);
    end
    if ~isempty(hd);
        hd.t = hd.tx - time_tx_min;
        hd.t = hd.t(hd.t>0,:);
        hd.pkg = hd.pkg(hd.t>0,:);
        hd.tx = hd.tx(hd.t>0,:);
        hd.rx = hd.rx(hd.t>0,:);
        hd.a_hd = hd.a_hd(hd.t>0,:);
        hd.a_eye = hd.a_eye(hd.t>0,:);
    end
    if ~isempty(tor);
        tor.t = tor.tx - time_tx_min;
        tor.t = tor.t(tor.t>0,:);
        tor.pkg = tor.pkg(tor.t>0,:);
        tor.tx = tor.tx(tor.t>0,:);
        tor.rx = tor.rx(tor.t>0,:);
        tor.a = tor.a(tor.t>0,:);
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
    if ~isempty(imgl);
        imgl.t = imgl.rx - time_rx_min;
        imgl.t = imgl.t(imgl.t>0,:);
        imgl.pkg = imgl.pkg(imgl.t>0,:);
        imgl.tx = imgl.tx(imgl.t>0,:);
        imgl.rx = imgl.rx(imgl.t>0,:);
        imgl.name = imgl.name(imgl.t>0,:);
    end
    if ~isempty(imgr);
        imgr.t = imgr.rx - time_rx_min;
        imgr.t = imgr.t(imgr.t>0,:);
        imgr.pkg = imgr.pkg(imgr.t>0,:);
        imgr.tx = imgr.tx(imgr.t>0,:);
        imgr.rx = imgr.rx(imgr.t>0,:);
        imgr.name = imgr.name(imgr.t>0,:);
    end
    if ~isempty(imu);
        imu.t = imu.rx - time_rx_min;
        imu.t = imu.t(imu.t>0,:);
        imu.pkg = imu.pkg(imu.t>0,:);
        imu.tx = imu.tx(imu.t>0,:);
        imu.rx = imu.rx(imu.t>0,:);
        imu.a = imu.a(imu.t>0,:);
        imu.f = imu.f(imu.t>0,:);
        imu.w = imu.w(imu.t>0,:);
        imu.B = imu.B(imu.t>0,:);
    end
    if ~isempty(hd);
        hd.t = hd.rx - time_rx_min;
        hd.t = hd.t(hd.t>0,:);
        hd.pkg = hd.pkg(hd.t>0,:);
        hd.tx = hd.tx(hd.t>0,:);
        hd.rx = hd.rx(hd.t>0,:);
        hd.a_hd = hd.a_hd(hd.t>0,:);
        hd.a_eye = hd.a_eye(hd.t>0,:);
    end
    if ~isempty(tor);
        tor.t = tor.rx - time_rx_min;
        tor.t = tor.t(tor.t>0,:);
        tor.pkg = tor.pkg(tor.t>0,:);
        tor.tx = tor.tx(tor.t>0,:);
        tor.rx = tor.rx(tor.t>0,:);
        tor.a = tor.a(tor.t>0,:);
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
