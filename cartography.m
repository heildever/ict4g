% ICT for geomatics - cartography exercises
clear all, close all, clc;
%%
% Ex.1 modulus of linear deformation
lambda(1) = 0;  % longtitude
phi = 30;   % latitude
i=1;
while lambda <= 3
    ml(i) = 0.9996*(1+(lambda(i)^2)*cos(phi).^2);
    i=i+1;
    lambda(i)= lambda(i-1)+0.5;    
end
plot(lambda(1:i-1),ml),title('modulus of linear deformation'), grid on;
%%
% Ex.2 Geo-2-Carto transformation using Hirvonen's Formula
% for P1 increase long by 6degrees to see 
clear all, clc ;
phi = [45,3,45.717];    % latitude
phi = degtorad(dms2degrees(phi));
lambda = [13,47,26.292];    % longtitude
lambda = degtorad(dms2degrees(lambda));
% ellipsoid parameters of WGS84
a = 6378137;
alfa = 1/298.257223563;
c = a*(1-alfa);
e = sqrt((a^2-c^2)/a^2);    % eccentricity of WGS84 ellipsoid
Rp = a^2/c;     % radius of polar curvature
??? lambda_mc = [12,27,8];   % central meridian of fuse 32
lambda_mc = degtorad(dms2degrees(lambda_mc));

??? e_prime = ????
??? lambda_prime = lambda-lambda_mc; 
v1 = sqrt(1+(e_prime^2)*cos(phi).^2);
Xi = atan(tan(phi)/cos(v1*lambda_prime));
v = sqrt(1+(e_prime^2)*cos(Xi).^2);

A1 = 1-(e^2)/4-3*(e^2)/64-5*(e^6)/256;
A2 = 3*(e^2)/4+3*(e^4)/32+45*(e^6)/1024;
A4 = 15*(e^4)/256+45*(e^6)/1024;
A6 = 35*(e^6)/3072;

x = Rp*asinh(cos(Xi)*tan(lambda_prime)/v) ;
y = a*(A1*Xi-A2*sin(2*Xi)+A4*sin(4*Xi)-A6*sin(6*Xi));
mc = 0.9996;    % contraction modulus 
East = x*mc+500;    % 500 = false east for fuse 32
North = y*mc;



