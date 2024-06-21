clear all; close all; clc;
%% Part 1A: Cylindrical Port (No Erosive Burning) CEA
% Uses O/F ratio of 4, metallization fraction of 25% Al, new density

% Fuel grain geometry
L_port = 35/100;                    % Length of propellant (m)
D_out = 6.6/100;                    % Outer propellant diameter (m)
D_port = 3/100;                     % Cylinder start diameter (m)    
rho_prop = 2828.3;                  % Propellant density (kg/m^3)

% Nozzle geometry
A_throat = 1.887/100^2;             % Nozzle throat area (m^2)
AR = 4;
A_exit = A_throat*AR;               % Nozzle exit area (m^2)
theta_exit = deg2rad(20);           % Nozzle exit angle (radians)

% Combustion gas properties (CEA Optimal)
gamma = 1.1565;
gamma_star = 1.1648;                  % Throat value
Mw = 24.56;                           % Molecular weight (kg/kg-mol)
Mw_star = 24.744;                     % Throat value  
T0 = 3133.3;                          % Flame temperature (K)

% Burn parameters
a = 0.132/100;                      % Burn multiplier (m/sec-kPa^n)
n = 0.16;                           % Burn exponent
M_crit = 0.3;                       % Critical Mach number
k = 0.2;                            % Mach scale factor

Ru = 8314.4126;                     % Universal gas constant (J/K-kg-mol)
Rg = Ru/Mw;                         % Specific gas constant (J/kg-K)
g0 = 9.8067;                        % Gravitational acceleration (m/s^2)

M_exit = mach_solver(AR,gamma,1,50);    % Mach number at nozzle exit

% Temperature at nozzle exit (K)
T_exit = T0 / (1 + ((gamma+1)/2)*M_exit^2);

% Velocity at nozzle exit (m/s)
V_exit = M_exit * sqrt(gamma * Rg * T_exit);

% Initial states
P0 = 101.325;                       % Chamber pressure (kPa)
r0 = D_port/2;                       % Port radius (m)
m0 = 0;                             % Propellant consumed (kg)

x0 = [P0; r0];
x_cylA(:,1) = x0;

dt = 0.001;                         % Time step

% Integration (Runge-Kutta 4)
i = 1;
while x_cylA(2) < D_out/2

    % Burn area (cm^2)
    A_burn = 2 * pi * x_cylA(2) * L_port;

    % Chamber volume (cm^3)
    Vc = pi * x_cylA(2)^2 * L_port;

    % K1
    [xdot] = cylinder1(x_cylA,a,n,rho_prop,Rg,T0,A_throat,gamma_star,A_burn,Vc,Mw,Mw_star);
    k_1 = xdot;
    xc_new = x_cylA + dt/2*k_1;

    % K2
    [xdot] = cylinder1(xc_new,a,n,rho_prop,Rg,T0,A_throat,gamma_star,A_burn,Vc,Mw,Mw_star);
    k_2 = xdot;
    xc_new = xc_new + dt/2*k_2;

    % K3
    [xdot] = cylinder1(xc_new,a,n,rho_prop,Rg,T0,A_throat,gamma_star,A_burn,Vc,Mw,Mw_star);
    k_3 = xdot;
    xc_new = xc_new + dt*k_3;

    % K4
    [xdot] = cylinder1(xc_new,a,n,rho_prop,Rg,T0,A_throat,gamma_star,A_burn,Vc,Mw,Mw_star);
    k_4 = xdot;

    % Chamber pressure
    x_cylA = x_cylA + (dt/6)*(k_1 + 2*k_2 + 2*k_3 + k_4);
    x_plot_cylA(i,:) = x_cylA;

    i = i+1;
end

% Time values
t_1A = zeros(1,length(x_plot_cylA));
for i = 1:length(t_1A)
    t_1A(i) = dt*i;
end

%% Part 2A: Cylindrical Port with Erosive Burning CEA
% The fuel grain geometry, nozzle geometry, and burn parameters are all the
% same as above for the erosive burn analysis.

% Initial states
x_eroA(:,1) = x0;

% Integration (Runge-Kutta 4)
i = 1;
while x_eroA(2) < D_out/2

    % Burn area (m^2)
    A_burn = 2 * pi * x_eroA(2) * L_port;

    % Chamber volume (m^3)
    Vc = pi * x_eroA(2)^2 * L_port;

    % K1
    [xdot] = erosive1(x_eroA,a,n,rho_prop,Rg,T0,A_throat,gamma_star,A_burn,Vc,k,M_crit,Mw,Mw_star);
    k_1 = xdot;
    xe_new = x_eroA + dt/2*k_1;

    % K2
    [xdot] = erosive1(xe_new,a,n,rho_prop,Rg,T0,A_throat,gamma_star,A_burn,Vc,k,M_crit,Mw,Mw_star);
    k_2 = xdot;
    xe_new = xe_new + dt/2*k_2;

    % K3
    [xdot] = erosive1(xe_new,a,n,rho_prop,Rg,T0,A_throat,gamma_star,A_burn,Vc,k,M_crit,Mw,Mw_star);
    k_3 = xdot;
    xe_new = xe_new + dt*k_3;

    % K4
    [xdot] = erosive1(xe_new,a,n,rho_prop,Rg,T0,A_throat,gamma_star,A_burn,Vc,k,M_crit,Mw,Mw_star);
    k_4 = xdot;

    % Chamber pressure and radius
    x_eroA = x_eroA + (dt/6)*(k_1 + 2*k_2 + 2*k_3 + k_4);
    x_plot_eroA(i,:) = x_eroA;

    i = i+1;
end

t21A = zeros(1,length(x_plot_eroA));
for i = 1:length(t21A)
    t21A(i) = dt*i;
end

%% Part 3A: Bates Grain CEA
% The fuel grain geometry, nozzle geometry, and burn parameters are all the
% same as above except for k.

% Mach scale factor for Bates grain
k1 = 0;

% Number of propellant sections
N_grain = 3;

% Initial states
x_batesA(:,1) = x0;

% Integration (Runge-Kutta 4)
i = 1;
while x_batesA(2) < D_out/2

    % Changing port diameter (m)
    D_port1 = x_batesA(2) * 2;

    % Radius of burnt propellant (m)
    s = x_batesA(2) - D_port/2;
    
    % Burn area (m^2)
    A_burn = N_grain * pi * ((D_out^2 - (D_port + 2*s)^2)/2 + (L_port/N_grain - 2*s)*(D_port + 2*s));
    
    % Chamber volume (m^3)
    Vc = (N_grain*pi)/4 * ((D_port + 2*s)^2 * (L_port/N_grain - 2*s) + (D_out^2 * 2*s));

    % K1
    [xdot] = bates1(x_batesA,a,n,rho_prop,Rg,T0,A_throat,gamma_star,A_burn,Vc,Mw,Mw_star);
    k_1 = xdot;
    xb_new = x_batesA + dt/2*k_1;

    % K2
    [xdot] = bates1(xb_new,a,n,rho_prop,Rg,T0,A_throat,gamma_star,A_burn,Vc,Mw,Mw_star);
    k_2 = xdot;
    xb_new = xb_new + dt/2*k_2;

    % K3
    [xdot] = bates1(xb_new,a,n,rho_prop,Rg,T0,A_throat,gamma_star,A_burn,Vc,Mw,Mw_star);
    k_3 = xdot;
    xb_new = xb_new + dt*k_3;

    % K4
    [xdot] = bates1(xb_new,a,n,rho_prop,Rg,T0,A_throat,gamma_star,A_burn,Vc,Mw,Mw_star);
    k_4 = xdot;

    % Chamber pressure and radius
    x_batesA = x_batesA + (dt/6)*(k_1 + 2*k_2 + 2*k_3 + k_4);
    x_plot_batesA(i,:) = x_batesA;

    i = i+1;
end

t31A = zeros(1,length(x_plot_batesA));
for i = 1:length(t31A)
    t31A(i) = dt*i;
end

%% Plot

% Chamber pressure profile
figure
plot(t_1A,x_plot_cylA(:,1),'LineWidth',1),hold on
plot(t21A,x_plot_eroA(:,1),'LineWidth',1)
plot(t31A,x_plot_batesA(:,1),'LineWidth',1)
% title('Chamber Pressure (kPa)')
xlabel('Time (s)')
ylabel('P_{0} (kPa)')
legend('Non-Erosive', 'Erosive', 'Bates-Grain','location', 'southeast')