function options = constraint_graph_optimise_set_options(options, field, val)

if nargin ~= 0 % change field
    assert(isfield(options, field), 'Invalid field')
    options.(field) = val;
else  % set default configuration
    options.iterations = 2;
    options.gateinnov = chi_square_bound(0.99, 6);
    options.gateresid = chi_square_bound(0.9, 6);
    options.gatetrust = chi_square_bound(0.9, 6);
    options.gateratio = 0.5;
    options.checkrank = true;
    options.verbose = 0; % verbose level: 0,1,2
end
