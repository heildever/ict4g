% ICT for geomatics - geodesy exercises
clear all, close all, clc;
%%
% Ex.1 ellipsoid parameter estimation
a = input('enter value of a');
alpha = input('enter value of alpha = ');
c = a*(1-alpha);
e = sqrt((a^2 - c^2)/a^2);
e_prime = sqrt((a^2 - c^2)/c^2);
disp(c),disp(e^2),disp(e_prime^2)

%%
% Ex.2 coordinate transformation from(phi,lambda,h) to (X, Y, Z)
% phi = input('enter phi value as an array = ');
% phi = degtorad(dms2degrees(phi));
% lambda = input('enter lambda value as an array = ');
% lambda = deg2rad(dms2degrees(lambda));
% h = input('enter h value as an array = ');
% W = sqrt(1-e^2*sin(phi));
% X = (a/W+h)*cos(phi)*cos(lambda);
% Y = (a/W+h)*cos(phi)*sin(lambda);
% Z = (a/W*(1-e^2)+h)*sin(phi);
% disp(X),disp(Y),disp(Z)

%%
% Ex.3 coordinate transformation from (X, Y, Z) to (phi,lambda,h)  
X = input('Enter X value = ');
Y = input ('Enter Y value = ');
Z = input ('Enter Z value = ');
lambda = atan(Y/X);
%step1
r = sqrt(X^2+Y^2) ;
phi = atan(Z/r);
%step2
W = sqrt(1-(e^2*sin(phi).^2));
%step3
N = a/W;
h = (X/(cos(phi)*cos(lambda)))-N; 
phi_prime = atan(Z/(r*(1-((e^2*N)/(N+h)))));
% first iteration 
W_2 = sqrt(1-(e^2*sin(phi_prime).^2));
N_2 = a/W_2;
h_2 = (X/(cos(phi)*cos(lambda)))-N_2;

phi_difference = phi-phi_prime;
h_difference = h-h_2;

    if  phi_difference <= 1e-8
        if h_difference <= 1e-8 
            W_2 = sqrt(1-e^2*sin(phi_prime));
            N_2 = a/W_2;
            h_2 = X/cos(phi_prime)*cos(lambda)-N_2; 
            phi_prime = atan(Z/r*(1-(e^2*N_2)/(N+h_2)));
        end
    end
    
%%    
clear all;  
% Ex.4 elevation and azimuth
% satellite coordinates
Xs = 15487292.829;
Ys = 6543538.932;
Zs = 20727274.429;

lat = [45,3,45.717];    % phi
lat = degtorad(dms2degrees(lat));
long = [13,47,26.292];    % lambda 
long = degtorad(dms2degrees(long));
% ellipsoid parameters WGS84
a = 6378137;
alpha = 1/298.257223563;
c = a*(1-alpha);
e = sqrt((a^2 - c^2)/a^2);
W = sqrt(1-(e^2*sin(lat).^2));
N = a/W;

Xi = N*cos(lat)*cos(long);
Yi = N*cos(lat)*sin(long);
Zi = N*(1-e^2)*sin(lat);
% geo coordinates of the station
R = [-sin(long) cos(long) 0;-sin(lat)*cos(long) -sin(lat)*sin(long) cos(lat); cos(lat)*cos(long) cos(lat)*sin(long) sin(lat)];
ENU = R*[Xs-Xi;Ys-Yi;Zs-Zi];
E = ENU(1);
N = ENU(2);
U = ENU(3);

elev = atan(U/sqrt(N^2+E^2));
azim = atan(E/N);
