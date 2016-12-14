function poses = sphereToEuclidean(azm, zen)
y = cos(zen);
x = sin(azm).*sin(zen);
z = -cos(azm).*sin(zen);
poses=[x y z];
return