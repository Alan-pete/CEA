function xdot = cylinder1(x,a,n,rho,Rg,T0,A_throat,gamma,A_burn,Vc,Mw,Mw_star)

P = x(1);

% Rate of change in chamber pressure (kPa/s)
Pdot = (A_burn * a * P^n)/Vc * (rho*Rg*T0/1000 - P) - (A_throat/Vc)*P*sqrt(gamma*Rg*T0*(Mw/Mw_star)*(2/(gamma+1))^((gamma+1)/(gamma-1)));

% Regression rate (cm/s)
rdot = a * P^n;

xdot = [Pdot; rdot];

end