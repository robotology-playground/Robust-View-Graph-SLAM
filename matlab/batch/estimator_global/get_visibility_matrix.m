function vis = get_visibility_matrix(impoints, siftind)

% vis : is the output visibility map

if nargin==2;
    
    x=[]; y=[]; feat=[];
    nimages = size(impoints.points,1);
    for i=1:nimages;
        x   =[x    impoints.index{i}]; % impoints.index{i} contains track numbers in the ith image
        y   =[y    repmat(i,size(impoints.index{i}))];
        feat=[feat siftind{i}];
    end
    vis = sparse(x,y,feat); 
    
else

    x=[]; y=[]; feat=[];
    nimages = size(impoints.points,1);
    for i=1:nimages;
        x   =[x    impoints.index{i}]; % impoints.index{i} contains track numbers in the ith image
        y   =[y    repmat(i,size(impoints.index{i}))];
    end
    vis = sparse(x,y,ones(size(y)));
    %vis = sparse(x,y,y);
    
end