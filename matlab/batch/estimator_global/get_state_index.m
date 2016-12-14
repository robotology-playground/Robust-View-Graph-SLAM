function idx=get_state_index(i)
% Get indices of (x,y,phi) for pose i in the state vector
idx=(1:6)+(i-1)*6;
