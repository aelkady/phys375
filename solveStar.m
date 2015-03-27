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
Kappa_H = (2.5e-32)*(Z/0.02)*((s(1)/(1e3))^(0.5))*(s(2)^(9)); 
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