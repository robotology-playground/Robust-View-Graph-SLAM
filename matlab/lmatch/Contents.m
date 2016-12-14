%LMATCH - MATCHING LINE SEGMENTS ACROSS MULTIPLE CALIBRATED IMAGES
%
%  - Coded by Tomas Werner with assistance of Andrew Zisserman, 2002.
%    Visual Geometry Group, Dept. of Engineering Science, University of Oxford, UK
%  - Improved by Tomas Werner, Center for Machine Perception, Prague, 2007.
%    werner@cmp.felk.cvut.cz
%
%  This is a public research software. Use for any other than research
%  purposes is prohibited. No responsibility for correct function 
%  of the software is accepted and no support is provided.
%-----------------------------------------------------------------------
%
%The MATLAB software 'LMATCH' allows matching line segments across
%multiple calibrated images. Its input are image bitmaps, camera
%projection matrices, and image line segments. Its output are line
%segment matches (tuples of image line segments, each from different
%image) and reconstructed 3D line segments.
%
%The matching is done in three steps:
%
% 1. Generating tentative matches.
%    A set of matches is generated, where each match is in accordance 
%    with the images (satisfies photommetric and geometric constraints),
%    but the matches may be inconsistent (in conflict) with each other.
%
%    This step is itself done in several substeps. In each substep,
%    matches are generated based on a given pair of views, so called
%    base view pair.
% 
% 2. Finding a consistent subset of matches.
%    A (sub)optimal consistent subset of tentative matches is found, 
%    and accepted as the final set of matches.
%
% 3. Reconstruction.
%    3D lines are reconstructed from the found lines segment matches.
%
%Compile MEXes and run 'lmatch_demo.m'.
%
%
%PUBLIC FUNCTIONS:
% lmatch_demo              - runs examples from ./examples/* and generates VRML file with 3D lines
% lmatch_detect_lines      - detects line segments in images (works only in Windows)
% lmatch_lineseg_orient    - orients line segments in a single image such that lighter part is always on one side
% lmatch_find_basepairs    - finds suitable base view pairs for generating tentative matches
% lmatch_options           - returns default options
% lmatch_generate          - generate tentative matches for a given base pair
% lmatch_resolve           - resolving tentative matches to final matches
% lmatch_reconstruct       - reconstructs 3D lines from line segment matches
% lmatch_vrml              - saves 3D lines to a vrml file
%
%See lmatch_memo.pdf for description of implementation details.
