function[function_rho_c, R_star, T_star, L_star, M_star, R, Rho, Temp, Mass, Lum] = getErrorInDensity(rho_c,T_c)
stefan_boltz = 5.670e-8; %sigma

%Fraction Things
X = 0.7; %1 - (2*10^-5);
X_CNO = 0.03*X;
Y = 0.28; %10^(-5);
Z = 1-X-Y; %10^(-5);

r_c = 1e-6;
R_sun = 6.955e8;
R_max = 50*R_sun;

M_c = ((4*pi) / 3)*(r_c^3)*(rho_c);

Epp_c = 1.07e-7*(rho_c / 1e5)*X^2*(T_c / 1e6)^4;
Ecno_c = 8.24e-26*(rho_c / 1e5)*X*X_CNO*(T_c / 1e6)^(19.9);
E_c = Epp_c + Ecno_c;

L_c = ((4*pi) / 3)*(r_c^3)*(rho_c)*E_c;

Kappa_es = 0.02*(1+X);
Kappa_ff = (1.0e24)*(1+X)*(Z+0.0001)*((rho_c/1e3)^(0.7))*(T_c^(-3.5));
Kappa_H = (2.5e-32)*(Z/0.02)*((rho_c/1e3)^(0.5)).*(T_c^(9));
Kappa = ((1/Kappa_H) + (1/max(Kappa_es, Kappa_ff))).^(-1);
Tau_c = Kappa*rho_c*r_c;

boundary_conditions_c = [rho_c, T_c, M_c, L_c, Tau_c];
options = odeset('Events', @starEvent);
[R, Star, ~, ~, ~, ~] = ode45(@solveStar, [r_c R_max], boundary_conditions_c, options);

Rho = Star(:,1);
Temp = Star(:,2);
Mass = Star(:,3);
Lum = Star(:,4);
Tau = Star(:,5);

[column,~] = size(R);
tauInfinity = Tau(column);
tauRstar = tauInfinity - (2/3);
[~,index] = min(abs(Tau-tauRstar));

R_star = R(index);
T_star = Temp(index);
L_star = Lum(index);
M_star = Mass(index);
function_rho_c = (L_star - 4*pi*stefan_boltz*((R_star)^2)*((T_star)^4))/((4*pi*stefan_boltz*((R_star)^2)*((T_star)^4)*L_star)^(1/2));
end