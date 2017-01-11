% a consule to compile and test mex files under construction

config_install;

% COMPILE

if (spinv)
	spinv_install();	
end

if (pwg)
    pwg_install();
end
