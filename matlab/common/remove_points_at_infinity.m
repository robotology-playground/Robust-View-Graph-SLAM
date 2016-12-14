function vis = remove_points_at_infinity(cpts1, cpts2, thres)

% matches should move at least pixel_disparity pixels
% the larget the pixel_disparity threshold is, the more stable the results
% and the estimates, but the less dense the result is and the shorter the
% observed range is.
if size(cpts1, 1) > 2
    if any(cpts1(3, :) ~= 1); % not in homo coordinates, it needs transpose
        cpts1 = cpts1';
        cpts2 = cpts2';
    end
end

diff1 = abs( cpts2(1, :) - cpts1(1, :) );
diff2 = abs( cpts2(2, :) - cpts1(2, :) );
vis = ( diff1 > thres ) | ( diff2 > thres );