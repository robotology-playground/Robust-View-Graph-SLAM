
function nrviews =  get_nrviews_in_each_track(impoints)

% counts the number of views in each track
% total number of tracks = impoints.pointnr
% indecies of tracks in the ith image = impoints.index{i}

nimages = size(impoints.points,1);
ntracks = impoints.pointnr;
nrviews = zeros(1, ntracks);
for i = 1:nimages;
    nrviews(impoints.index{i}) = nrviews(impoints.index{i}) + 1;
end
