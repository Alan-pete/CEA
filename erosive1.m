function xdot = erosive1(x,a,n,rho,Rg,T0,A_throat,gamma,A_burn,Vc,k,M_crit,Mw,Mw_star)

P = x(1);
r = x(2);

% Port area
A_port = pi * r^2;

% Port Mach number
M_port = mach_solver(A_port/A_throat, gamma, 0, 50);

% Regression rate (m/s)
rdot = ((1 + k*(M_port/M_crit)) * a*P^n) / (1 + k);

% Rate of change in chamber pressure (kPa/s)
Pdot = A_burn*(((1 + k*(M_port/M_crit)) * a*P^n) / (1 + k))*(rho*Rg*T0/1000 - P)/Vc - (A_throat/Vc)*P*sqrt(gamma*Rg*T0*(Mw/Mw_star)*(2/(gamma+1))^((gamma+1)/(gamma-1)));

xdot = [Pdot; rdot];

% Propellant massflow (kg/s)
%mdot_prop = rho * A_burn * rdot;

end