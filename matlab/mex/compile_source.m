%--------------------------------------------------------------------------
function compile_source(d, include, output, mex_src, libs)
fprintf(['Compiling ' output '\n']);
s = sprintf ('mex %s -DDLONG -O %s -output %s %s', d, include, output, mex_src);
s = [s ' ' libs];
eval(s);