function result = debug_conditional(level)
%function result = debug_conditional(level)
%
% INPUT: 
%   level - set level of conditional, three levels available
%
% OUTPUT:
%   result - 1 if condition accepted, 0 if condition failed
%
% Permit adjustable levels of debug testing, so that expensive tests may be
% turned off when faster operation is required. To set a test-level, call
% the function without a return argument as 
%
%       debug_conditional(level)
% 
% The level is one of four values
%       0 - permit no tests
%       1 - permit cheap tests only 
%       2 - permit non-expensive tests (cheap, moderate)
%       3 - permit all tests (cheap, moderate, expensive)
%
% To use the conditional, call the function without an argument
%
%       if debug_conditional(level)
%           do_test
%       end
%
% In this mode of operation, the level values correspond to 
%       1 - cheap test
%       2 - moderate test
%       3 - expensive test
%
% Tim Bailey 2011.
persistent LEVEL
if level < 0 || level > 3
    error('Invalid test level, must be in {0,1,2,3}')
end

if nargout == 0
    LEVEL = level;
else
    if level == 0
        error('Invalid test level, tests must be in {1,2,3}')
    end
    if isempty(LEVEL)
        warning('Debug level has not been set. Defaulting to level = 1')
        LEVEL = 1;
    end
    result = level <= LEVEL;
end
