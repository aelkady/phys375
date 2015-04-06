format long

G = 6.67384e-11; %gravitational constant
h_bar = 1.054571726e-34; %reduced planck constant
m_e = 9.10938291e-31; %mass of an electron
m_p = 1.67262178e-27; %mass of a proton
stefan_boltz = 5.670e-8; %sigma
c = 2.99792458e8; %speed of light
a = (4*stefan_boltz) / c; %a thing
k = 1.3806488e-23; %boltzmann constant
Gamma = (5/3);

X = 0.7; %1 - (2*10^-5);
X_CNO = 0.03*X;
Y = 0.28;
Z = 1-X-Y; 
mu = (2*X + 0.75*Y + 0.5*Z)^(-1);

%T_c = 35e7;
M_sun = 1.989e30;
L_sun = 3.846e26;
R_sun = 6.95800e8;

eps_abs = 1e-5;
eps_step = 1e-5;

A = 0.2e6:5e5:10e6;
B = 10e6:2.5e5:100e6;
T_c = cat(2,A,B);
length = size(T_c, 2);
Luminosity = ones(1, length);
T_star = ones(1, length);

for j=1:length 
    rho_c_min = 300;
    rho_c_max = 500000;
    function_rho_c_new = 100000;
    i=0;
    while ((abs(real(function_rho_c_new)) > eps_abs) && (i<200))
        rho_c_new = real((rho_c_min + rho_c_max)/2);
        [function_rho_c_min, ~, ~, ~, ~, ~, ~, ~, ~, ~] = getErrorInDensity(rho_c_min,T_c(1, j));
        [function_rho_c_max, ~, ~, ~, ~, ~, ~, ~, ~, ~] = getErrorInDensity(rho_c_max,T_c(1, j));
        [function_rho_c_new, R_star_new, T_star_new, L_star_new, M_star_new, R, Rho, Temp, Mass, Lum] = getErrorInDensity(rho_c_new,T_c(1, j));
%         if(real(function_rho_c_min)*real(function_rho_c_max) > 0)
%            rho_c_new = 0;
%            R_star_new = 0;
%            T_star_new = 0;
%            L_star_new = 0;
%            M_star_new = 0;
%            R = 0;
%            Rho = 0;
%            Temp = 0;
%            Mass = 0;
%            Lum = 0;
%            break;
%         end
        if (real(function_rho_c_new) == 0)
           break;
        elseif ( function_rho_c_new > 0 )
           rho_c_max = rho_c_new;
        else
           rho_c_min = rho_c_new;
        end
        if ( rho_c_max - rho_c_min < eps_step )
            if ( abs( real(function_rho_c_min) ) < abs( real(function_rho_c_max) ) && abs( real(function_rho_c_min) ) < eps_abs )
                rho_c_new = rho_c_min;
                break;
            elseif ( abs( real(function_rho_c_max) ) < eps_abs )
                rho_c_new = rho_c_max;
                break;
            end
        end
        i=i+1;
    end
    i = 0;
    if((L_star_new == 0) || (abs(real(function_rho_c_new)) > 7.5))
        Luminosity(1, j) = nan;
    else   
        Luminosity(1, j) = real(L_star_new/L_sun);
    end;
    if((T_star_new == 0) || (abs(real(function_rho_c_new)) > 7.5))
        T_star(1, j) = nan;
    else   
        T_star(1, j) = real(T_star_new);
    end;
end

plot(log10(T_star), log10(Luminosity), '*b');
set(gca,'xdir','reverse');
xlim([3 4.15])
ylim([-6 6])
xlim([3 4.1])
ylim([-6 5])
title('Main Sequence of Stars')
xlabel('Log Base 10 of Temperature')
ylabel('Log Base 10 of Luminosity Divided by the Luminosity of the Sun')

