
function CMD;
global delta_t;
%CONSTANTS
G = 6.67384*10^(-11); %gravitational constant
h_bar = 1.054571726*10^(-34); %reduced planck constant
m_e = 9.10938291*10^(-31); %mass of an electron
m_p = 1.67262178*10^(-27); %mass of a proton
stefan_boltz = 5.670*10^(-8); %sigma
c = 299792458; %speed of light
a = (4*stefan_boltz) / c; %a thing
k = 1.3806488*10^(-23); %boltzmann constant
Gamma = (5/3);

%Fraction Things
X = 0.73; %1 - (2*10^-5);
X_CNO = 0.03*X;
Y = 0.25; %10^(-5);
Z = 0.02; %10^(-5);

mu = (2*X + 0.75*Y + 0.5*Z)^(-1);

r_c = 0.001;
r_max_1 = 10^7;
r_max_2 = 10^9;
%r = [r_c, 10^10];
r = linspace(r_c, 10^9, 10);
T_c = 15*10^6;
rho_c = 1.622*10^5;
M_c = ((4*pi) / 3)*(r_c^3)*(rho_c);

Epp_c = 1.07*10^(-7)*(rho_c / 10^5)*X^2*(T_c / 10^6)^4;
Ecno_c = 8.24*10^(-26)*(rho_c / 10^5)*X*X_CNO*(T_c / 10^6)^(-19.9);
E_c = Epp_c + Ecno_c;

L_c = ((4*pi) / 3)*(r_c^3)*(rho_c)*E_c;

boundary_conditions_c = [rho_c, T_c, M_c, L_c, 10^5];

%options = odeset('RelTol',1e-4,'AbsTol',[1e-5 1e-5 1e-5 1e-5 1e-5]);
%for rMax = [r_c:r_max]
for rMax = [10^7:10^9]
    options = odeset('Events', @starEvent);
    [R, Star, TE, YE, IE, sol] = ode45(@solveStar, [r_c rMax], boundary_conditions_c, options);
    if(isEmpty(IE))
        fprintf('sup');
    else
        fprintf('found it');
        break
    end
end
%Star(:,1) = rho, Star(:,2) = T
Kappa_es = ones(10,1)*(0.02*(1+X));
Kappa_ff = (1.0.*10^24).*(Z+0.0001).*(Star(:,1).^(0.7)).*(Star(:,2).^(-3.5));
Kappa_H = (2.5.*10.^(-32)).*(Z/0.02).*(Star(:,1).^(0.5)).*(Star(:,2).^(9));
Kappa = ((1./Kappa_H) + (1./max(Kappa_es, Kappa_ff))).^(-1);
Star(:,1)
%hold all;
%plot(R, log10(Kappa_es));
%plot(R, log10(Kappa_ff));
%plot(R, log10(Kappa_H));
%plot(R, log10(Kappa));
%hold off;

%plot(R, Star(:,1));
%title('Rho');
%figure();
%plot(R, Star(:,2));
%title('Temp');
%figure();
%plot(R, Star(:,3));
%title('Mass');
%figure();
%plot(R, Star(:,4));
%title('Lum');
%figure();
%plot(R, Star(:,5));
%title('Tau');
%plot(R, Kappa);
%title('Kappa');

function ds = solveStar(r, s)
global delta_t;
%SHHHHH Constants are here again
G = 6.67384*10^(-11); %gravitational constant
h_bar = 1.054571726*10^(-34); %reduced planck constant
m_e = 9.10938291*10^(-31); %mass of an electron
m_p = 1.67262178*10^(-27); %mass of a proton
stefan_boltz = 5.670*10^(-8); %sigma
c = 299792458; %speed of light
a = (4*stefan_boltz) / c; %a thing
k = 1.3806488*10^(-23); %boltzmann constant
Gamma = (5/3);

%Fraction Things
X = 1 - (2*10^-5);
X_CNO = 0.03*X;
Y = 10^(-5);
Z = 10^(-5);
mu = (2*X + 0.75*Y + 0.5*Z)^(-1);

%FUNCTION PROPER START

%ORDER: rho:1, T:2, M:3, L:4, Tau:5
ds = zeros(5, 1);

%HELPER EQUATIONS
Kappa_es = 0.02*(1+X);
Kappa_ff = (1.0*10^24)*(Z+0.0001)*(s(1)^(0.7))*(s(2)^(-3.5));
Kappa_H = (2.5*10^(-32))*(Z/0.02)*(s(1)^(0.5))*(s(2)^(9)); 
Kappa = ((1/Kappa_H) + (1/max(Kappa_es, Kappa_ff)))^(-1);

%drho/dr
dP_dT = (s(1)*k) / (mu*m_p) + (4/3)*a*(s(2)^3);
dP_drho = (((3*pi^2)^(2/3))/3)*((h_bar^2)/m_e)*((s(1)/m_p)^(2/3)) + (k*s(2)/(mu*m_p));
ds(1) = -1*(((G*s(3)*s(1)) / (r^2)) + dP_dT*ds(2)) / dP_drho;

%dT/dr ... has the min stuff
Pressure = (((3*pi^2)^(2/3)) / 5)*((h_bar^2) / m_e) * (s(1)/m_p)^(5/3) + s(1)*(k*s(2) / mu*m_p) + (1/3)*a*(s(2)^4);
dT_L = (3*Kappa*s(1)*s(4)) / (16*pi*a*c*(s(2)^3)*(r^2));
dT_R = (1 - (1/Gamma))*(s(2) / Pressure)*((G*s(3)*s(1)) / r^2);
ds(2) = -1*min(dT_L, dT_R);

%dM/dr 
ds(3) = 4*pi*r^2*s(1);

%dL/dr
Epp = 1.07*10^(-7)*(s(1) / 10^5)*X^2*(s(2) / 10^6)^4;
Ecno = 8.24*10^(-26)*(s(1) / 10^5)*X*X_CNO*(s(2) / 10^6)^(-19.9);
E = Epp + Ecno;
ds(4) = 4*pi*r^2*s(1)*E;

%dTau/dr
ds(5) = Kappa*s(1);
delta_t = (Kappa*(s(1)^2)) / (abs(ds(1)));

function[value, isTerminal, direction] = starEvent(r, s)
global delta_t;
thresh = 0.001;
value = s(5);
isTerminal = (delta_t < thresh);
direction = 0;

