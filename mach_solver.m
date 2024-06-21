function [M] = mach_solver(AR,gamma,big,n)
% A,A_star
% A/A_star = AR just changed for final exam
% determining the first guess for Mach number
if big == 1
    m_init = 3.;
else
    m_init = 0.01;
end

% starting a list for the iterations of the Mach number, starts with the
% initial guess
M_n = zeros(1, n);
M_n(1) = m_init;

error = 1;

% for loop to iterate through
while error > 0.001
    for i = 2:n
        % numerator in the Newton solver
        F_m = 1/M_n(i-1) * ((2/(gamma + 1)) * (1 + (gamma - 1)/2 * M_n(i-1)^2))^((gamma + 1)/(2*(gamma - 1))) - AR; %A/A_star;
    
        % denominator in the Newton solver
        dF_dm1 = 2^((1 - 3*gamma)/(2 - 2*gamma)) * (M_n(i-1)^2 - 1)/((M_n(i-1)^2) * (2 + (M_n(i-1)^2)*(gamma - 1)));
        dF_dm2 = ((1 + (gamma - 1)/2 * (M_n(i-1)^2)) / (gamma + 1))^((gamma + 1)/(2*(gamma - 1)));
    
        % full Newton solver
        M_n(i) = M_n(i-1) - (F_m / (dF_dm1 * dF_dm2));
    
        % error to determine when to drop from loop
        error = abs(1/M_n(i) * ((2/(gamma + 1)) * (1 + (gamma - 1)/2 * M_n(i)^2))^((gamma + 1)/(2*(gamma - 1))) - AR) / AR; %A/A_star) / (A/A_star);

    end
end

M = M_n(n);





