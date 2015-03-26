
function CMD;
global delta_t;
%CONSTANTS
G = 6.67384e-11; %gravitational constant
h_bar = 1.054571726e-34; %reduced planck constant
m_e = 9.10938291e-31; %mass of an electron
m_p = 1.67262178e-27; %mass of a proton
stefan_boltz = 5.670e-8; %sigma
c = 299792458; %speed of light
a = (4*stefan_boltz) / c; %a thing
k = 1.3806488e-23; %boltzmann constant
Gamma = (5/3);


%Fraction Things
X = 0.7; %1 - (2*10^-5);
X_CNO = 0.03*X;
Y = 0.28; %10^(-5);
Z = 1-X-Y; %10^(-5);

mu = (2*X + 0.75*Y + 0.5*Z)^(-1);

r_c = 1e-6;
R_sun = 6.955e8;
R_max = 50*R_sun;
T_c = 8.23e6;
t = cputime;
i = 1;
rho_c = (300:100:500000);
sizeTest = size(rho_c);
rho_c_matrix = ones(sizeTest(2),1);

for rho_c = (300:100:500000)
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
    [R, Star, TE, YE, IE, sol] = ode45(@solveStar, [r_c R_max], boundary_conditions_c, options);

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
    rho_c_matrix(i,1) = function_rho_c;
    if(i ~= 1)
        if (rho_c_matrix(i,1) * rho_c_matrix(i-1,1) < 0)
            break;
        end
    end
    i = i + 1;
end

T_star
rho_c
T_c
R_star
M_star
L_star

plot(R, Rho);
title('Rho');
figure();
plot(R, Temp);
title('Temp');
figure();
plot(R, Mass);
title('Mass');
figure();
plot(R, Lum);
title('Lum');
figure();
plot(R, Tau);
title('Tau');
figure();

% [cols_size row_size] = size(R);
% Kappa_es = ones(cols_size,1)*(0.02*(1+X));
% Kappa_ff = (1.0.*10^24).*(1+X).*(Z+0.0001).*((Star(:,1)./1e3).^(0.7)).*(Star(:,2).^(-3.5));
% Kappa_H = (2.5.*10.^(-32)).*(Z/0.02).*((Star(:,1)./1e3).^(0.5)).*(Star(:,2).^(9));
% Kappa = ((1./Kappa_H) + (1./max(Kappa_es, Kappa_ff))).^(-1);
% figure();
% hold all;
% plot(R, log10(Kappa_es));
% plot(R, log10(Kappa_ff));
% plot(R, log10(Kappa_H));
% plot(R, log10(Kappa));
% hold off;
% title('Kappa');

function ds = solveStar(r, s)
global delta_t;
%SHHHHH Constants are here again
G = 6.67384e-11; %gravitational constant
h_bar = 1.054571726e-34; %reduced planck constant
m_e = 9.10938291e-31; %mass of an electron
m_p = 1.67262178e-27; %mass of a proton
stefan_boltz = 5.670e-8; %sigma
c = 2.99792458e8; %speed of light
a = (4*stefan_boltz) / c; %a thing
k = 1.3806488e-23; %boltzmann constant
Gamma = (5/3);

%Fraction Things
X = 0.7; %1 - (2*10^-5);
X_CNO = 0.03*X;
Y = 0.28;
Z = 1-X-Y; 
mu = (2*X + 0.75*Y + 0.5*Z)^(-1);

%FUNCTION PROPER START

%ORDER: rho:1, T:2, M:3, L:4, Tau:5
ds = zeros(5, 1);

%HELPER EQUATIONS
Kappa_es = 0.02*(1+X);
Kappa_ff = (1.0e24)*(1+X)*(Z+0.0001)*((s(1)/1e3)^(0.7))*(s(2)^(-3.5));
Kappa_H = (2.5e-32)*(Z/0.02)*((s(1)/1e3)^(0.5))*(s(2)^(9)); 
Kappa = ((1/Kappa_H) + (1/max(Kappa_es, Kappa_ff)))^(-1);


%dT/dr ... has the min stuff
Pressure = ((((3*(pi^2))^(2/3)) / 5)*((h_bar^2) / m_e) * (s(1)/m_p)^(5/3)) + ((s(1)*k*s(2)) / (mu*m_p)) + (1/3)*a*(s(2)^4);
dT_L = (3*Kappa*s(1)*s(4)) / (16*pi*a*c*(s(2)^3)*(r^2));
dT_R = (1 - (1/Gamma))*(s(2) / Pressure)*((G*s(3)*s(1)) / r^2);
ds(2) = -1*min(abs(dT_L), abs(dT_R));



%drho/dr
dP_dT = ((s(1)*k) / (mu*m_p)) + (4/3)*a*(s(2)^3);
dP_drho = (((3*pi^2)^(2/3))/3)*((h_bar^2)/(m_e*m_p))*((s(1)/m_p)^(2/3)) + (k*s(2)/(mu*m_p));
ds(1) = -1*(((G*s(3)*s(1)) / (r^2)) + dP_dT*ds(2)) / dP_drho;



%dM/dr 
ds(3) = 4*pi*r^2*s(1);

%dL/dr
Epp = (1.07e-7)*(s(1) / 1e5)*X^2*(((s(2) / 1e6))^4);
Ecno = (8.24e-26)*(s(1) / 1e5)*X*X_CNO*((s(2) / 1e6)^(19.9));
E = Epp + Ecno;
ds(4) = 4*pi*(r^2)*s(1)*E;

%dTau/dr
ds(5) = Kappa*s(1);
delta_t = (Kappa*(s(1)^2)) / (abs(ds(1)));

function[value, isTerminal, direction] = starEvent(~, s)
global delta_t;
thresh = 1e-3;
if((delta_t <= thresh || (s(3) > 1.989e33)))
    value = 0;
else
    value = 1;
end;
isTerminal = 1;
direction = 0;

