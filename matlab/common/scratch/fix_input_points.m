

% fix inputs into 3XN form
% Tariq Abuhashim - August 2014, iCub

function [x1,x2]=fix_input_points(x1,x2)


% get inputs size
[m1,n1]=size(x1);
[m2,n2]=size(x2);

% check taht both have the same dimensions
if m1~=m2 || n1~=n2;
    error(['input points have different dimensions: ',...
        'x1 is ',num2str(m1),'x',num2str(n1),' and ',...
        'x2 is ',num2str(m2),'x',num2str(n2)]);
end    

% check if at least one dimension is 2 or 3
if m1==2 || m1==3 || n1==2 || n1==3;
    % check if at least one dimension is larger than 3
    if m1>3 || n1>3;
        % check if transpose needed
        if m1>n1;
            x1=x1'; x2=x2';
        end
        % check if not homogenous
        [m,~]=size(x1);
        if m==2;
            x1=pextend(x1);
            x2=pextend(x2);
        end
    else
        error(['not enough input points: ',...
            'x1 is ',num2str(m1),'x',num2str(n1),' and ',...
            'x2 is ',num2str(m2),'x',num2str(n2)]);
    end
else
    error(['input points have wrong dimensions: ',...
        'x1 is ',num2str(m1),'x',num2str(n1),' and ',...
        'x2 is ',num2str(m2),'x',num2str(n2),'. ',...
        'rows or cols should have dimensions 2 or 3']);
end

% check if last row is ones
if sum(x1(3,:)==1)~=size(x1,2);
    error(['this function only works for the case of 2D points, ',...
        'last row should have all ones']);
end