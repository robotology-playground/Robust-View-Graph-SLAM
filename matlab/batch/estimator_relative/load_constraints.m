function [C, scan, kpts, ncam, npts] = load_constraints(options)

disp('  ') ;
load(strcat(options.save, '/constraints')) ;
n = 0 ;
for i = 1 : size(scan, 2) ;
    if ~isempty(scan{i}) ;
        n = n + 1 ;
    end ;
end ;
disp(['loaded ' num2str(n) ' x scan']) ;
n = 0 ;
for i = 1 : size(kpts, 2) ;
    if ~isempty(kpts{i}) ;
        n = n + 1 ;
    end ;
end ;
disp(['loaded ' num2str(n) ' x kpts']) ;
n = length(C) ;
disp(['loaded ' num2str(n) ' x constraints']) ;
disp('Press any key to continue, or Ctrl+c to exit ...') ;

pause;