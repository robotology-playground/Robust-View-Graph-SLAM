function str = zero_padding(number, strlength)

% pads the array of numbers "number" by zeros such that all elements have
% the same length.
%
% Tariq Abuhashim

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