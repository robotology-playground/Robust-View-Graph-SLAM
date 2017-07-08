
function out = read_calibration_file(FileName)
%out = read_calibration_file(FileName)
% Reads the calibration file "FileName" of HRP2 stereo with large baseline. Toulouse Calibration.
%
% Tariq Abuhashim, 2017.
%
% Auckland, New Zealand

fid = fopen(FileName, 'r');
if fid < 0, error('Cannot open file'); end

out = [];
len = 30; % locate memory first
part = zeros(1, len);
ipart = 0;
while 1  % Infinite loop
	s = fgets(fid);
	if ischar(s)
		data = sscanf(s, '%f', 1);
    	if length(data) == 1
      		ipart = ipart + 1;
      		part(:, ipart) = data;
      		if ipart == len
        		out = cat(2, out, part);
        		ipart = 0;
      		end
    	end
  	else  % End of file:
    	break;
  	end
end
out = cat(2, out, part(:, 1:ipart));
