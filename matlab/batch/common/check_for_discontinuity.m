function check_for_discontinuity(v, S)
if v'*(S\v) > chi_square_bound(.999, length(v)) % check for large discrepancy
    warning('NUMERICAL:discontinuity', 'Possible discontinuity found')
end
