%% Part 1: Cylindrical Port (No Erosive Burning)
close all; clear all; clc;

% Fuel grain geometry
L_port = 35/100;                    % Length of propellant (m)
D_out = 6.6/100;                    % Outer propellant diameter (m)
D_port = 3/100;                     % Cylinder start diameter (m)    
rho_prop = 1260;                    % Propellant density (kg/m^3)

% Nozzle geometry
A_throat = 1.887/100^2;             % Nozzle throat area (m^2)
AR = 4;
A_exit = A_throat*AR;               % Nozzle exit area (m^2)
theta_exit = deg2rad(20);           % Nozzle exit angle (radians)

% Combustion gas properties
gamma = 1.18;
Mw = 23;                            % Molecular weight (kg/kg-mol)
T0 = 2900;                          % Flame temperature (K)

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
T_exit = T0 / (1 + ((gamma-1)/2)*M_exit^2);

% Velocity at nozzle exit (m/s)
V_exit = M_exit * sqrt(gamma * Rg * T_exit);

% Initial states
P0 = 101.325;                       % Chamber pressure (kPa)
r0 = D_port/2;                       % Port radius (m)
m0 = 0;                             % Propellant consumed (kg)

x0 = [P0; r0];
x_cyl(:,1) = x0;
mdot_cyl(1) = 0;
rdot_cyl(1) = a*P0^n;
mdot_choke_cyl(1) = A_throat * sqrt((gamma/Rg) * (2/(gamma+1))^((gamma+1)/(gamma-1))) * (P0/sqrt(T0));
m_cyl(1) = 0;
p_cyl(1) = P0;
thrust_cyl(1) = 0;
impulse_cyl(1) = 0;
Isp_cyl(1) = 0;

dt = 0.001;                         % Time step

% Integration (Runge-Kutta 4)
i = 1;
while x_cyl(2) < D_out/2

    % Burn area (cm^2)
    A_burn = 2 * pi * x_cyl(2) * L_port;

    % Chamber volume (cm^3)
    Vc = pi * x_cyl(2)^2 * L_port;

    % K1
    [xdot] = cylinder(x_cyl,a,n,rho_prop,Rg,T0,A_throat,gamma,A_burn,Vc);
    k_1 = xdot;
    xc_new = x_cyl + dt/2*k_1;

    % K2
    [xdot] = cylinder(xc_new,a,n,rho_prop,Rg,T0,A_throat,gamma,A_burn,Vc);
    k_2 = xdot;
    xc_new = xc_new + dt/2*k_2;

    % K3
    [xdot] = cylinder(xc_new,a,n,rho_prop,Rg,T0,A_throat,gamma,A_burn,Vc);
    k_3 = xdot;
    xc_new = xc_new + dt*k_3;

    % K4
    [xdot] = cylinder(xc_new,a,n,rho_prop,Rg,T0,A_throat,gamma,A_burn,Vc);
    k_4 = xdot;

    % Chamber pressure
    x_cyl = x_cyl + (dt/6)*(k_1 + 2*k_2 + 2*k_3 + k_4);
    x_plot_cyl(i,:) = x_cyl;
    
    % Propellant massflow (kg/s)
    mdot_cyl(i+1) = rho_prop * A_burn * (a*x_cyl(1)^n);

    % Regression rate (m/s)
    rdot_cyl(i+1) = a*x_cyl(1)^n;

    % Choking massflow (kg/s)
    mdot_choke_cyl(i+1) = A_throat * sqrt((gamma/Rg) * (2/(gamma+1))^((gamma+1)/(gamma-1))) * (x_cyl(1)/sqrt(T0));

    % Mass depletion (kg)
    m_cyl(i+1) = m_cyl(i) + mdot_cyl(i+1)*dt;

    % Exit Pressure (kPa)
    p_cyl(i+1) = x_cyl(1) / (1 + ((gamma-1)/2)*M_exit^2)^(gamma/(gamma-1));

    % Thrust (N)
    thrust_cyl(i+1) = mdot_cyl(i+1)*V_exit + (p_cyl(i+1) - P0)*A_exit;

    % Total impulse (Ns)
    impulse_cyl(i+1) = impulse_cyl(i) + thrust_cyl(i+1)*dt;

    % Specific Impulse (sec)
    Isp_cyl(i+1) = thrust_cyl(i+1) / (mdot_cyl(i+1)*g0);

    i = i+1;
end

Isp_cyl_avg = mean(Isp_cyl);

% Time values
t_1 = zeros(1,length(x_plot_cyl));
for i = 1:length(t_1)
    t_1(i) = dt*i;
end

t = zeros(1,length(x_plot_cyl)+1);
for i = 1:length(t)
    t(i) = dt*i;
end

%% Part 2: Cylindrical Port with Erosive Burning
% The fuel grain geometry, nozzle geometry, and burn parameters are all the
% same as above for the erosive burn analysis.

% Initial states
x_ero(:,1) = x0;
mdot_ero(1) = 0;
M_port_ero(1) = mach_solver((pi*r0^2)/A_throat,gamma,0.1,50);
rdot_ero(1) = ((1 + k*(M_port_ero(1)/M_crit)) * a*P0^n) / (1 + k);
mdot_choke_ero(1) = A_throat * sqrt((gamma/Rg) * (2/(gamma+1))^((gamma+1)/(gamma-1))) * (P0/sqrt(T0));
m_ero(1) = 0;
p_ero(1) = P0;
thrust_ero(1) = 0;
impulse_ero(1) = 0;
Isp_ero(1) = 0;

% Integration (Runge-Kutta 4)
i = 1;
while x_ero(2) < D_out/2

    % Burn area (m^2)
    A_burn = 2 * pi * x_ero(2) * L_port;

    % Chamber volume (m^3)
    Vc = pi * x_ero(2)^2 * L_port;

    % K1
    [xdot] = erosive(x_ero,a,n,rho_prop,Rg,T0,A_throat,gamma,A_burn,Vc,k,M_crit);
    k_1 = xdot;
    xe_new = x_ero + dt/2*k_1;

    % K2
    [xdot] = erosive(xe_new,a,n,rho_prop,Rg,T0,A_throat,gamma,A_burn,Vc,k,M_crit);
    k_2 = xdot;
    xe_new = xe_new + dt/2*k_2;

    % K3
    [xdot] = erosive(xe_new,a,n,rho_prop,Rg,T0,A_throat,gamma,A_burn,Vc,k,M_crit);
    k_3 = xdot;
    xe_new = xe_new + dt*k_3;

    % K4
    [xdot] = erosive(xe_new,a,n,rho_prop,Rg,T0,A_throat,gamma,A_burn,Vc,k,M_crit);
    k_4 = xdot;

    % Chamber pressure and radius
    x_ero = x_ero + (dt/6)*(k_1 + 2*k_2 + 2*k_3 + k_4);
    x_plot_ero(i,:) = x_ero;
    
    % Port Mach number
    M_port_ero(i+1) = mach_solver((pi*x_ero(2)^2)/A_throat,gamma,0.1,50);

    % Regression rate (m/s)
    rdot_ero(i+1) = ((1 + k*(M_port_ero(i+1)/M_crit)) * a*x_ero(1)^n) / (1 + k);

    % Propellant massflow (kg/s)
    mdot_ero(i+1) = rho_prop * A_burn * rdot_ero(i+1);

    % Choking massflow (kg/s)
    mdot_choke_ero(i+1) = A_throat * sqrt((gamma/Rg) * (2/(gamma+1))^((gamma+1)/(gamma-1))) * (x_ero(1)/sqrt(T0));

    % Mass depletion (kg)
    m_ero(i+1) = m_ero(i) + mdot_ero(i+1)*dt;

    % Exit Pressure (kPa)
    p_ero(i+1) = x_ero(1) / (1 + ((gamma-1)/2)*M_exit^2)^(gamma/(gamma-1));

    % Thrust (N)
    thrust_ero(i+1) = mdot_ero(i+1)*V_exit + (p_ero(i+1) - P0)*A_exit;

    % Total impulse (Ns)
    impulse_ero(i+1) = impulse_ero(i) + thrust_ero(i+1)*dt;

    % Specific Impulse (sec)
    Isp_ero(i+1) = thrust_ero(i+1) / (mdot_ero(i+1)*g0);

    i = i+1;
end

Isp_ero_avg = mean(Isp_ero);

t21 = zeros(1,length(x_plot_ero));
for i = 1:length(t21)
    t21(i) = dt*i;
end

t2 = zeros(1,length(x_plot_ero)+1);
for i = 1:length(t2)
    t2(i) = dt*i;
end

%% Part 3: Bates Grain
% The fuel grain geometry, nozzle geometry, and burn parameters are all the
% same as above except for k.

% Mach scale factor for Bates grain
k1 = 0;

% Number of propellant sections
N_grain = 3;

% Initial states
x_bates(:,1) = x0;
mdot_bates(1) = 0;
rdot_bates(1) = a*P0^n;
mdot_choke_bates(1) = A_throat * sqrt((gamma/Rg) * (2/(gamma+1))^((gamma+1)/(gamma-1))) * (P0/sqrt(T0));
m_bates(1) = 0;
p_bates(1) = P0;
thrust_bates(1) = 0;
impulse_bates(1) = 0;
Isp_bates(1) = 0;

% Integration (Runge-Kutta 4)
i = 1;
while x_bates(2) < D_out/2

    % Changing port diameter (m)
    D_port1 = x_bates(2) * 2;

    % Radius of burnt propellant (m)
    s = x_bates(2) - D_port/2;
    
    % Burn area (m^2)
    A_burn = N_grain * pi * ((D_out^2 - (D_port + 2*s)^2)/2 + (L_port/N_grain - 2*s)*(D_port + 2*s));
    
    % Chamber volume (m^3)
    Vc = (N_grain*pi)/4 * ((D_port + 2*s)^2 * (L_port/N_grain - 2*s) + (D_out^2 * 2*s));

    % K1
    [xdot] = bates(x_bates,a,n,rho_prop,Rg,T0,A_throat,gamma,A_burn,Vc);
    k_1 = xdot;
    xb_new = x_bates + dt/2*k_1;

    % K2
    [xdot] = bates(xb_new,a,n,rho_prop,Rg,T0,A_throat,gamma,A_burn,Vc);
    k_2 = xdot;
    xb_new = xb_new + dt/2*k_2;

    % K3
    [xdot] = bates(xb_new,a,n,rho_prop,Rg,T0,A_throat,gamma,A_burn,Vc);
    k_3 = xdot;
    xb_new = xb_new + dt*k_3;

    % K4
    [xdot] = bates(xb_new,a,n,rho_prop,Rg,T0,A_throat,gamma,A_burn,Vc);
    k_4 = xdot;

    % Chamber pressure and radius
    x_bates = x_bates + (dt/6)*(k_1 + 2*k_2 + 2*k_3 + k_4);
    x_plot_bates(i,:) = x_bates;
    
    % Propellant massflow (kg/s)
    mdot_bates(i+1) = rho_prop * A_burn * (a*x_bates(1)^n);

    % Regression rate (m/s)
    rdot_bates(i+1) = a*x_bates(1)^n;

    % Choking massflow (kg/s)
    mdot_choke_bates(i+1) = A_throat * sqrt((gamma/Rg) * (2/(gamma+1))^((gamma+1)/(gamma-1))) * (x_bates(1)/sqrt(T0));

    % Mass depletion (kg)
    m_bates(i+1) = m_bates(i) + mdot_bates(i+1)*dt;

    % Exit Pressure (kPa)
    p_bates(i+1) = x_bates(1) / (1 + ((gamma-1)/2)*M_exit^2)^(gamma/(gamma-1));

    % Thrust (N)
    thrust_bates(i+1) = mdot_bates(i+1)*V_exit + (p_bates(i+1) - P0)*A_exit;

    % Total impulse (Ns)
    impulse_bates(i+1) = impulse_bates(i) + thrust_bates(i+1)*dt;

    % Specific Impulse (sec)
    Isp_bates(i+1) = thrust_bates(i+1) / (mdot_bates(i+1)*g0);

    i = i+1;
end

Isp_bates_avg = mean(Isp_bates);

t31 = zeros(1,length(x_plot_bates));
for i = 1:length(t31)
    t31(i) = dt*i;
end

t3 = zeros(1,length(x_plot_bates)+1);
for i = 1:length(t3)
    t3(i) = dt*i;
end

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
gamma = 1.156;
Mw = 24.547;                            % Molecular weight (kg/kg-mol)
T0 = 3140.9;                          % Flame temperature (K)

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
mdot_cylA(1) = 0;
rdot_cylA(1) = a*P0^n;
mdot_choke_cylA(1) = A_throat * sqrt((gamma/Rg) * (2/(gamma+1))^((gamma+1)/(gamma-1))) * (P0/sqrt(T0));
m_cylA(1) = 0;
p_cylA(1) = P0;
thrust_cylA(1) = 0;
impulse_cylA(1) = 0;
Isp_cylA(1) = 0;

dt = 0.001;                         % Time step

% Integration (Runge-Kutta 4)
i = 1;
while x_cylA(2) < D_out/2

    % Burn area (cm^2)
    A_burn = 2 * pi * x_cylA(2) * L_port;

    % Chamber volume (cm^3)
    Vc = pi * x_cylA(2)^2 * L_port;

    % K1
    [xdot] = cylinder(x_cylA,a,n,rho_prop,Rg,T0,A_throat,gamma,A_burn,Vc);
    k_1 = xdot;
    xc_new = x_cylA + dt/2*k_1;

    % K2
    [xdot] = cylinder(xc_new,a,n,rho_prop,Rg,T0,A_throat,gamma,A_burn,Vc);
    k_2 = xdot;
    xc_new = xc_new + dt/2*k_2;

    % K3
    [xdot] = cylinder(xc_new,a,n,rho_prop,Rg,T0,A_throat,gamma,A_burn,Vc);
    k_3 = xdot;
    xc_new = xc_new + dt*k_3;

    % K4
    [xdot] = cylinder(xc_new,a,n,rho_prop,Rg,T0,A_throat,gamma,A_burn,Vc);
    k_4 = xdot;

    % Chamber pressure
    x_cylA = x_cylA + (dt/6)*(k_1 + 2*k_2 + 2*k_3 + k_4);
    x_plot_cylA(i,:) = x_cylA;
    
    % Propellant massflow (kg/s)
    mdot_cylA(i+1) = rho_prop * A_burn * (a*x_cylA(1)^n);

    % Regression rate (m/s)
    rdot_cylA(i+1) = a*x_cylA(1)^n;

    % Choking massflow (kg/s)
    mdot_choke_cylA(i+1) = A_throat * sqrt((gamma/Rg) * (2/(gamma+1))^((gamma+1)/(gamma-1))) * (x_cylA(1)/sqrt(T0));

    % Mass depletion (kg)
    m_cylA(i+1) = m_cylA(i) + mdot_cylA(i+1)*dt;

    % Exit Pressure (kPa)
    p_cylA(i+1) = x_cylA(1) / (1 + ((gamma-1)/2)*M_exit^2)^(gamma/(gamma-1));

    % Thrust (N)
    thrust_cylA(i+1) = mdot_cylA(i+1)*V_exit + (p_cylA(i+1) - P0)*A_exit;

    % Total impulse (Ns)
    impulse_cylA(i+1) = impulse_cylA(i) + thrust_cylA(i+1)*dt;

    % Specific Impulse (sec)
    Isp_cylA(i+1) = thrust_cylA(i+1) / (mdot_cylA(i+1)*g0);

    i = i+1;
end

Isp_cyl_avgA = mean(Isp_cylA);

% Time values
t_1A = zeros(1,length(x_plot_cylA));
for i = 1:length(t_1A)
    t_1A(i) = dt*i;
end

tA = zeros(1,length(x_plot_cylA)+1);
for i = 1:length(tA)
    tA(i) = dt*i;
end

%% Part 2A: Cylindrical Port with Erosive Burning CEA
% The fuel grain geometry, nozzle geometry, and burn parameters are all the
% same as above for the erosive burn analysis.

% Initial states
x_eroA(:,1) = x0;
mdot_eroA(1) = 0;
M_port_eroA(1) = mach_solver((pi*r0^2)/A_throat,gamma,0.1,50);
rdot_eroA(1) = ((1 + k*(M_port_eroA(1)/M_crit)) * a*P0^n) / (1 + k);
mdot_choke_eroA(1) = A_throat * sqrt((gamma/Rg) * (2/(gamma+1))^((gamma+1)/(gamma-1))) * (P0/sqrt(T0));
m_eroA(1) = 0;
p_eroA(1) = P0;
thrust_eroA(1) = 0;
impulse_eroA(1) = 0;
Isp_eroA(1) = 0;

% Integration (Runge-Kutta 4)
i = 1;
while x_eroA(2) < D_out/2

    % Burn area (m^2)
    A_burn = 2 * pi * x_eroA(2) * L_port;

    % Chamber volume (m^3)
    Vc = pi * x_eroA(2)^2 * L_port;

    % K1
    [xdot] = erosive(x_eroA,a,n,rho_prop,Rg,T0,A_throat,gamma,A_burn,Vc,k,M_crit);
    k_1 = xdot;
    xe_new = x_eroA + dt/2*k_1;

    % K2
    [xdot] = erosive(xe_new,a,n,rho_prop,Rg,T0,A_throat,gamma,A_burn,Vc,k,M_crit);
    k_2 = xdot;
    xe_new = xe_new + dt/2*k_2;

    % K3
    [xdot] = erosive(xe_new,a,n,rho_prop,Rg,T0,A_throat,gamma,A_burn,Vc,k,M_crit);
    k_3 = xdot;
    xe_new = xe_new + dt*k_3;

    % K4
    [xdot] = erosive(xe_new,a,n,rho_prop,Rg,T0,A_throat,gamma,A_burn,Vc,k,M_crit);
    k_4 = xdot;

    % Chamber pressure and radius
    x_eroA = x_eroA + (dt/6)*(k_1 + 2*k_2 + 2*k_3 + k_4);
    x_plot_eroA(i,:) = x_eroA;

    % Port Mach number
    M_port_eroA(i+1) = mach_solver((pi*x_eroA(2)^2)/A_throat,gamma,0.1,50);

    % Regression rate (m/s)
    rdot_eroA(i+1) = ((1 + k*(M_port_eroA(i+1)/M_crit)) * a*x_eroA(1)^n) / (1 + k);
    
    % Propellant massflow (kg/s)
    mdot_eroA(i+1) = rho_prop * A_burn * rdot_eroA(i+1);

    % Choking massflow (kg/s)
    mdot_choke_eroA(i+1) = A_throat * sqrt((gamma/Rg) * (2/(gamma+1))^((gamma+1)/(gamma-1))) * (x_eroA(1)/sqrt(T0));

    % Mass depletion (kg)
    m_eroA(i+1) = m_eroA(i) + mdot_eroA(i+1)*dt;

    % Exit Pressure (kPa)
    p_eroA(i+1) = x_eroA(1) / (1 + ((gamma-1)/2)*M_exit^2)^(gamma/(gamma-1));

    % Thrust (N)
    thrust_eroA(i+1) = mdot_eroA(i+1)*V_exit + (p_eroA(i+1) - P0)*A_exit;

    % Total impulse (Ns)
    impulse_eroA(i+1) = impulse_eroA(i) + thrust_eroA(i+1)*dt;

    % Specific Impulse (sec)
    Isp_eroA(i+1) = thrust_eroA(i+1) / (mdot_eroA(i+1)*g0);

    i = i+1;
end

Isp_ero_avgA = mean(Isp_eroA);

t21A = zeros(1,length(x_plot_eroA));
for i = 1:length(t21A)
    t21A(i) = dt*i;
end

t2A = zeros(1,length(x_plot_eroA)+1);
for i = 1:length(t2A)
    t2A(i) = dt*i;
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
mdot_batesA(1) = 0;
rdot_batesA(1) = a*P0^n;
mdot_choke_batesA(1) = A_throat * sqrt((gamma/Rg) * (2/(gamma+1))^((gamma+1)/(gamma-1))) * (P0/sqrt(T0));
m_batesA(1) = 0;
p_batesA(1) = P0;
thrust_batesA(1) = 0;
impulse_batesA(1) = 0;
Isp_batesA(1) = 0;

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
    [xdot] = bates(x_batesA,a,n,rho_prop,Rg,T0,A_throat,gamma,A_burn,Vc);
    k_1 = xdot;
    xb_new = x_batesA + dt/2*k_1;

    % K2
    [xdot] = bates(xb_new,a,n,rho_prop,Rg,T0,A_throat,gamma,A_burn,Vc);
    k_2 = xdot;
    xb_new = xb_new + dt/2*k_2;

    % K3
    [xdot] = bates(xb_new,a,n,rho_prop,Rg,T0,A_throat,gamma,A_burn,Vc);
    k_3 = xdot;
    xb_new = xb_new + dt*k_3;

    % K4
    [xdot] = bates(xb_new,a,n,rho_prop,Rg,T0,A_throat,gamma,A_burn,Vc);
    k_4 = xdot;

    % Chamber pressure and radius
    x_batesA = x_batesA + (dt/6)*(k_1 + 2*k_2 + 2*k_3 + k_4);
    x_plot_batesA(i,:) = x_batesA;
    
    % Propellant massflow (kg/s)
    mdot_batesA(i+1) = rho_prop * A_burn * (a*x_batesA(1)^n);

    % Regression rate (m/s)
    rdot_batesA(i+1) = a*x_batesA(1)^n;

    % Choking massflow (kg/s)
    mdot_choke_batesA(i+1) = A_throat * sqrt((gamma/Rg) * (2/(gamma+1))^((gamma+1)/(gamma-1))) * (x_batesA(1)/sqrt(T0));

    % Mass depletion (kg)
    m_batesA(i+1) = m_batesA(i) + mdot_batesA(i+1)*dt;

    % Exit Pressure (kPa)
    p_batesA(i+1) = x_batesA(1) / (1 + ((gamma-1)/2)*M_exit^2)^(gamma/(gamma-1));

    % Thrust (N)
    thrust_batesA(i+1) = mdot_batesA(i+1)*V_exit + (p_batesA(i+1) - P0)*A_exit;

    % Total impulse (Ns)
    impulse_batesA(i+1) = impulse_batesA(i) + thrust_batesA(i+1)*dt;

    % Specific Impulse (sec)
    Isp_batesA(i+1) = thrust_batesA(i+1) / (mdot_batesA(i+1)*g0);

    i = i+1;
end

Isp_bates_avgA = mean(Isp_batesA);

t31A = zeros(1,length(x_plot_batesA));
for i = 1:length(t31A)
    t31A(i) = dt*i;
end

t3A = zeros(1,length(x_plot_batesA)+1);
for i = 1:length(t3A)
    t3A(i) = dt*i;
end

%% Plots

% Chamber pressure profile
figure
plot(t_1,x_plot_cyl(:,1),'LineWidth',1), hold on
plot(t21,x_plot_ero(:,1),'LineWidth',1)
plot(t31,x_plot_bates(:,1),'LineWidth',1)
plot(t_1A,x_plot_cylA(:,1),'LineWidth',1)
plot(t21A,x_plot_eroA(:,1),'LineWidth',1)
plot(t31A,x_plot_batesA(:,1),'LineWidth',1)
% title('Chamber Pressure (kPa)')
xlabel('Time (s)')
ylabel('P_{0} (kPa)')
legend('Non-Erosive', 'Erosive', 'Bates-Grain','Non-Erosive CEA', 'Erosive CEA', 'Bates-Grain CEA','location', 'southeast')

% Regression rate profile
figure
plot(t,rdot_cyl,'LineWidth',1), hold on
plot(t2,rdot_ero,'LineWidth',1)
plot(t3,rdot_bates,'LineWidth',1)
plot(tA,rdot_cylA,'LineWidth',1)
plot(t2A,rdot_eroA,'LineWidth',1)
plot(t3A,rdot_batesA,'LineWidth',1)
% title('Regression Rate')
xlabel('Time (s)')
ylabel('rdot (m/s)')
legend('Non-Erosive', 'Erosive', 'Bates-Grain','Non-Erosive CEA', 'Erosive CEA', 'Bates-Grain CEA','location', 'southeast')

% Thrust profile
figure
plot(t,thrust_cyl,'LineWidth',1), hold on
plot(t2,thrust_ero,'LineWidth',1)
plot(t3,thrust_bates,'LineWidth',1)
plot(tA,thrust_cylA,'LineWidth',1)
plot(t2A,thrust_eroA,'LineWidth',1)
plot(t3A,thrust_batesA,'LineWidth',1)
% title('Thrust')
xlabel('Time (s)')
ylabel('Thrust (N)')
legend('Non-Erosive', 'Erosive', 'Bates-Grain','Non-Erosive CEA', 'Erosive CEA', 'Bates-Grain CEA','location', 'southeast')

% Total impulse profile
figure
plot(t,impulse_cyl,'LineWidth',1), hold on
plot(t2,impulse_ero,'LineWidth',1)
plot(t3,impulse_bates,'LineWidth',1)
plot(tA,impulse_cylA,'LineWidth',1)
plot(t2A,impulse_eroA,'LineWidth',1)
plot(t3A,impulse_batesA,'LineWidth',1)
% title('Total Impulse')
xlabel('Time (s)')
ylabel('Impulse (Ns)')
legend('Non-Erosive', 'Erosive', 'Bates-Grain','Non-Erosive CEA', 'Erosive CEA', 'Bates-Grain CEA','location', 'northwest')

% Specific Impulse
figure
plot(t,Isp_cyl,'LineWidth',1), hold on
plot(t2,Isp_ero,'LineWidth',1)
plot(t3,Isp_bates,'LineWidth',1)
plot(tA,Isp_cylA,'LineWidth',1)
plot(t2A,Isp_eroA,'LineWidth',1)
plot(t3A,Isp_batesA,'LineWidth',1)
% title('Specific Impulse')
xlabel('Time (s)')
ylabel('I_{sp} (sec)')
legend('Non-Erosive', 'Erosive', 'Bates-Grain','Non-Erosive CEA', 'Erosive CEA', 'Bates-Grain CEA','location', 'southeast')

