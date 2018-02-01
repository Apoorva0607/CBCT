% phantom_def.m
% Defines 3D Shepp-Logan phantom paramters.
% Note: the phantom is defined in units of cm for lengths

function E=sheppdef3D()

% define parameters for the 10 objects comprising the 3D S-L phantom
% xhat = x-coordinate of center 
% yhat = y-coordinate of center 
% zhat = z-coordinate of center
% a = first half axis length
% b = second half axis length
% c = third half axis length
% theta  = angle between n3 and the z-axis (in degrees)
% phi = angle between n3 and the x-axis (in degrees)
% psi = angle between n1 and e1 -- see picture in course notes (in degrees)
% mu = x-ray attentuation coefficient 

xhat = [0 0 0.22 -0.22 0 0 0 -0.08 0 0.06]'*10;
yhat = [0 -0.0184 0 0 0.35 0.1 -0.1 -0.605 -0.605 -0.605]'*10;
zhat = [0 0 0 0 0 0 0 0 0 0]'*10;
a = [0.69 0.6624 0.11 0.16 0.21 0.046 0.046 0.046 0.023 0.023]'*10;
b = [0.92 0.874 0.31 0.41 0.25 0.046 0.046 0.023 0.023 0.046]'*10;
c = [.9 .88 .21 .22 .35 .046 .02 .02 .1 .1]'*10;
theta = zeros(10,1);
phi = [0 0 -18 18 0 0 0 0 0 0]';  
psi = zeros(10,1);
mu =  [2 -0.98 -0.02 -0.02 -0.01 0.01 0.01 0.01 0.01 0.01]';

E = [xhat,yhat,zhat,a,b,c,theta,phi,psi,mu];


