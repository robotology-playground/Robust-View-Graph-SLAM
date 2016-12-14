%lmatch_options  Returns default options.
%
% opt = lmatch_options returns default options for functions 'lmatch_generate' and 'lmatch_resolve'.
% Any option can be changed, e.g. 'opt.Calibration=3'.
%
% Options affecting function 'lmatch_generate' only:
% - 'DispRange' ... displacement search range along epipolars (default: Inf pixels).
%      Says how far two matching line segments can be from each other. This constraint is imposed not 
%      on each view pair from the collection, but only to:
%       - base view pairs
%       - pairs of nearest (in terms of view distance D) view pair during searching other views
%      Lower (but realistic) value implies lower runtime and better results.
% - 'MaxViewDistance' ... the maximal view distance between query view (in which a match has no segment) and 
%     the nearest (in terms of view distance D) view where the match has a segment (default inf).
%
% Options affecting function 'lmatch_generate', and function 'lmatch_resolve' with 'Merging=1':
% - 'Calibration' ... calibration level: 1 = (oriented) projective, 2 = affine, 3 = metric (default 3).
%      Higher level often implies better matching results.
% - 'NCCWindow' ... vector of parameters of correlation window used to compute similarity of image line segments:
%   - NCCWindow(1) = half side of the correlation window perpendicular to the line segment (default: 6 pixels)
%   - NCCWindow(2) = half side of the correlation window parallel to the line segment (default: 7 pixels)
%   - NCCWindow(3) = distance of the centers of neighboring windows on the line segment (default: 3 pixels)
%   - NCCWindow(4) = 0 means window centered on the line; 1 means window centered on the line and also on its both sides (default: 1)
% - 'NCCThreshold' ... vector of thresholds on NCC to accept the match:
%   - NCCThreshold(1) = threshold on NCC value (default: 0.6)
%   - NCCThreshold(2) = threshold on number of segment's pixels having NCC higher than NCCThresh(1) (default: 10 pixels)
% - 'ReprojResidual' ... 2-vector, reprojection residual threshold for linear/nonlinear estimation (default: [5 1.5] pixels).
%      Matching line segments must satisfy reprojection test up to this accuracy. Increase for inaccurate calibration.
%
% Options affecting function 'lmatch_resolve' only:
% 'Ordering' ... 0 = no ordering constraint, 1 = ord. constraint (slow but accurate), 2 = ord. constraint fast but less accurate (default: 0)
%    Often (i.e., if no thin objects are in foreground), line segments obey ordering constraint.
%    Imposing eliminates some mismatches but does not match possible thin objects in the foreground.
%    There are two algorithms for imposing ordering constraints, one slower and one faster.
%    Resolving with Ordering~=0 is quadratic in number of line segments, hence time consuming for large image sets.
% 'Merging' ... 0 = merging not applied, 1 = merging applied (default: 0).
%    Due to errors in line detection, image line segments are often fragmented, and this fragmentation is in general
%    different in each view. Some of these this fragmented segments can be merged again based on evidence in 3D.
%    This process can be time consuming for large number of views.
% 'MergeResidual' ... scalar, max. residual for merging fragmented line segments (default: 1)

function opt = lmatch_options

opt.DispRange        = Inf;
opt.MaxViewDistance  = Inf;

opt.Calibration      = 1;
opt.ReprojResidual   = [5 1.5];
opt.NCCWindow        = [6 7 3 1];
opt.NCCThreshold     = [.6 10];
 
opt.Ordering         = logical(0);
opt.Merging          = logical(0);
opt.MergeResidual    = 1;

return