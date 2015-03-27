
format long

T_c = 8.23e6;
M_sun = 1.989e30;
L_sun = 3.846e26;
R_sun = 6.95800e8;

eps_abs = 1e-2;
eps_step = 1e-2;
rho_c_min = 300;
rho_c_max = 500000;
function_rho_c_min = 300;
function_rho_c_max = 500000;

while (rho_c_max - rho_c_min >= eps_step || ( abs( function_rho_c_min ) >= eps_abs && abs( function_rho_c_max )  >= eps_abs ) )
    rho_c_new = (rho_c_min + rho_c_max)/2;
    [function_rho_c_min, ~, ~, ~, ~, ~, ~, ~, ~, ~] = getErrorInDensity(rho_c_min,T_c);
    [function_rho_c_max, ~, ~, ~, ~, ~, ~, ~, ~, ~] = getErrorInDensity(rho_c_max,T_c);
    [fucntion_rho_c_new, R_star_new, T_star_new, L_star_new, M_star_new, R, Rho, Temp, Mass, Lum] = getErrorInDensity(rho_c_new,T_c);
    if (fucntion_rho_c_new == 0)
       break;
    elseif ( function_rho_c_min*fucntion_rho_c_new < 0 )
       rho_c_max = rho_c_new;
    else
       rho_c_min = rho_c_new;
    end
end

%e = t - cputime
T_star_new
rho_c_new
T_c
R_star_new/R_sun
M_star_new/M_sun
L_star_new/L_sun

hold on
plot(R./R_star_new, Rho./rho_c_new, '-k');
title('Rho');
plot(R./R_star_new, Temp./T_c, '-r');
title('Temp');
plot(R./R_star_new, Mass./M_star_new,'-g');
title('Mass');
plot(R./R_star_new, Lum./L_star_new,'-b');
title('Lum');

[cols_size,~] = size(R);
Kappa_es = ones(cols_size,1)*(0.02*(1+X));
Kappa_ff = (1.0.*10^24).*(1+X).*(Z+0.0001).*((Rho./1e3).^(0.7)).*(Temp.^(-3.5));
Kappa_H = (2.5.*10.^(-32)).*(Z/0.02).*((Rho./1e3).^(0.5)).*(Temp.^(9));
Kappa = ((1./Kappa_H) + (1./max(Kappa_es, Kappa_ff))).^(-1);
figure();
hold all;
plot(R, log10(Kappa_es));
plot(R, log10(Kappa_ff));
plot(R, log10(Kappa_H));
plot(R, log10(Kappa));
hold off;
title('Kappa');

