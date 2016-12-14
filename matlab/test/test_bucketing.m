function buckets = test_bucketing(pts, dims, img)

bucket_width = dims(1);
bucket_height = dims(2);

if size(pts,1) > 2; pts = pts'; end;

% find max values
u_max = max(pts(1,:));
v_max = max(pts(2,:));

% allocate number of buckets needed
bucket_cols = floor(u_max/bucket_width ) + 1;
bucket_rows = floor(v_max/bucket_height) + 1;
  
% assign matches to their buckets
u = floor(pts(1,:)/bucket_width ) + 1; u_ind = unique(u);
%if any(u_ind<0)
%    warning('Negative backets have been removed.');
%    u_ind = u_ind(u_ind>0);
%end
v = floor(pts(2,:)/bucket_height) + 1; v_ind = unique(v);
%if any(v_ind<0)
%    warning('Negative backets have been removed.');
%    v_ind = v_ind(v_ind>0);
%end

% fill the buckets
%buckets=cell(length(u_ind)*length(v_ind),1);
buckets = cell(bucket_cols*bucket_rows, 1);
for i = 1:length(u_ind);%bucket_cols  % u
    for j = 1:length(v_ind);%bucket_rows  % v
        buckets{bucket_rows*(u_ind(i)-1)+v_ind(j)} = find(u==u_ind(i) & v==v_ind(j));
    end
end

% plot the buckets
if nargin > 2
    imshow(img);hold on;
    for i = 1:bucket_cols
        plot([bucket_width*i bucket_width*i],[1 size(img,1)],'linewidth',2);
    end
    for i = 1:bucket_rows
        plot([1 size(img,2)],[bucket_height*i bucket_height*i],'linewidth',2);
    end
    for i = 1:size(buckets,1)
        %plot(pts(1,buckets{i}),pts(2,buckets{i}),'+','color',[rand rand rand],...
        %    'markersize',2,'linewidth',2);
        plot(pts(1,buckets{i}),pts(2,buckets{i}),'+','markersize',5,'linewidth',2);
    end
    drawnow;
end