function invK = inverse_K(K)

% inverse of calibration  matrix (K)
invK = zeros(size(K));

tau = K(1,1)/K(2,2);
invK(1,1) = 1/(tau*K(2,2));
invK(1,2) = 0;
invK(1,3) = -K(1,3)/(tau*K(2,2));
invK(2,1) = 0;
invK(2,2) = 1/K(2,2);
invK(2,3) = -K(2,3)/K(2,2);
invK(3,1) = 0;
invK(3,2) = 0;
invK(3,3) = 1;