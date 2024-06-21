close all; clear all; clc;

D_t = 1.55;                 % Throat diameter (cm)
D_e = 3.1;                  % Exit diameter (cm)
R1 = 0.75 * D_t;            % Radius of curvature of nozzle contraction (cm)
L_n = 2.13;                 % Nozzle length (cm)
theta = deg2rad(20);        % Nozzle angle (rad)

gamma = 1.18;               % Ratio of specific heats
A_t = pi * D_t^2/4;         % Nozzle throat area (cm^2)
A_e = pi * D_e^2/4;         % Nozzle exit area (cm^2)

% Area ratio
AR = A_e/A_t;

% Mach at exit
M_e = mach_solver(AR,gamma,1,50);

% Angle at exit
nu_e_rad = sqrt((gamma+1)/(gamma-1)) * atan(sqrt((gamma-1)/(gamma+1) * (M_e^2 - 1))) - atan(sqrt(M_e^2 - 1));
nu_e = rad2deg(nu_e_rad);

% Theta max
theta_max = nu_e/2;

% New minimum length
L_new1 = 1/2*(sqrt(AR) - 1)*D_t + R1*(1/cosd(theta_max) - 1);
L_new2 = tand(theta_max);

L_new = L_new1 / L_new2

% Factor of safety
FOS = (nu_e/2 - 20) / 20
