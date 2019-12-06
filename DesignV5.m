%% ASEN3113 Design Lab
% Tristan Liu, Colin Ruark, Evan Welch, Maria Callas, Roland Bailey

% fat refresh button
clear; close all; clc

%% general

q_BL = [88 63 75.5 11]; % IR backload at WS, SS, equinox, eclipse [W/m^2]
opPower = 20; % operational power [W]

T_op = 25 + 273; % operational temperature [K]
T_max = 30 + 273; % max instrument temperature [K]
T_su = -40 + 273; % survival temperature [K]

eps = 0.85; % emissivity
alpha = 0.2; % absorptivity 
sigma = 5.67e-8; % Stefan-Boltzmann constant [W/(m^2 K^4)]

tilt = deg2rad(23.5); % Earth's tilt at solstice (rad)

q_eq = 1361; % flux at equinox
q_w = 1384; % flux at winter solstice
q_s = 1338; % flux at summer solstice
q_ec = 0; % flux at eclipse

q_flux = [q_w q_s q_eq q_ec]; % solar flux at WS, SS, equinox, eclipse [W/m^2]

%% area calculations

% designed at 30 deg C such that surface area is smallest it can be

for i=1:4
    A_mat(i) = opPower / (eps*sigma*(T_max)^4 - alpha*q_flux(i) - eps*q_BL(i));
end

% calculated required areas at...
A_w = A_mat(1); % WS
A_s = A_mat(2); % SS
A_eq = A_mat(3); % equinox
A_ec = A_mat(4); % eclipse

% The highest possible values occur across the board when the satellite is
% in equinox, which makes sense because that is when the solar radiation is
% directly perpendicular to the surface area of the radiator. 

%% time in eclipse (at equinox)

R_E = 6378; % radius of Earth [km]
h_sat = 42164; % altitude of a GEO orbit [km]
R_sat = R_E + h_sat; % radius of GEO orbit, aka location of satellite [km]
d_E = 2* R_E; % diameter of Earth [km]

% the length of the orbit in which the satellite is eclipsed by the Earth
% is equal to the diameter of Earth (*assumption*)

% arc length = radius * theta
thetaE = d_E / R_sat; % duration the satelite is eclipsed [rad]

thetaFrac = thetaE / (2*pi); % fraction of the orbit that sat is eclipsed [-]

tEclipsed = thetaFrac * 24; % duration the satellite is eclipsed by the Earth at equinox [hr]

%% operational heater power as a function of time

% designed to be held at 25 deg C

periodE = 0:0.1:24; % period of one Earth day [hr]
theta = 0:(2*pi)/240:(2*pi); % radian location of satellite [rad]
A_rot_eq = A_eq .* sin(theta); % radiator area with incidence angle

% winter solstice
for i=1:length(theta)
    if (theta(i) <= pi) % in sunlight
        Q_inf_w_op(i) = A_eq * (eps*sigma*T_op^4  - eps*q_BL(1)) - A_rot_eq(i)*alpha*q_w - 20;
        T_w_op(i) = ((20 + A_rot_eq(i)*alpha*q_w + A_eq*eps*q_BL(1)) / (eps*sigma*A_eq)).^(1/4) - 273;
        Q_solar_rad_w(i) = sin(theta(i))*q_w;
    elseif (theta(i) > pi)
        Q_inf_w_op(i) = A_eq * (eps*sigma*T_op^4 - eps*q_BL(1)) - 20;
        Q_solar_rad_w(i) = 0;
        T_w_op(i) = ((20 + A_eq*eps*q_BL(1)) / (eps*sigma*A_eq)).^(1/4) - 273;
    end
end

figure(1)
grid on
yyaxis left
plot(periodE,Q_inf_w_op,'LineWidth',3)
xlabel('Time (Hr)')
ylabel('Power (W)')
yyaxis right
plot(periodE,T_w_op,'LineWidth',3)
xlabel('Time (Hr)')
ylabel('Temperature (^{\circ}C)')
title('Operational Winter Solstice')

% summer solstice
for i=1:length(theta)
    if (theta(i) <= pi) % in sunlight
        Q_inf_s_op(i) = A_eq * (eps*sigma*T_op^4  - eps*q_BL(2)) - A_rot_eq(i)*alpha*q_s - 20;
        T_s_op(i) = ((20 + A_rot_eq(i)*alpha*q_s + A_eq*eps*q_BL(2)) / (eps*sigma*A_eq)).^(1/4) - 273;
        Q_solar_rad_s(i) = sin(theta(i))*q_s;
    elseif (theta(i) > pi)
        Q_inf_s_op(i) = A_eq * (eps*sigma*T_op^4 - eps*q_BL(2)) - 20;
        Q_solar_rad_s(i) = 0;
        T_s_op(i) = ((20 + A_eq*eps*q_BL(2)) / (eps*sigma*A_eq)).^(1/4) - 273;
    end
end

figure(2)
grid on
yyaxis left
plot(periodE,Q_inf_s_op,'LineWidth',3)
xlabel('Time (Hr)')
ylabel('Power (W)')
yyaxis right
plot(periodE,T_s_op,'LineWidth',3)
xlabel('Time (Hr)')
ylabel('Temperature (^{\circ}C)')
title('Operational Summer Solstice')

% equinox
for i=1:length(theta)
    if (theta(i) > (pi - thetaE/2)) && (theta(i) < (pi + thetaE/2)) % eclipse
        Q_inf_e_op(i) = A_eq * (sigma*eps*T_op^4 - eps*q_BL(4)) - 20;
        T_e_op(i) = ((A_eq*eps*q_BL(4) + 20) / (sigma*eps*A_eq)).^(1/4) - 273;
        Q_solar_rad_eq(i) = 0;
    elseif (theta(i) <= pi) % in sunlight
        Q_inf_e_op(i) = A_eq * (eps*sigma*T_op^4  - eps*q_BL(3)) - A_rot_eq(i)*alpha*q_eq - 20;
        T_e_op(i) = ((20 + A_rot_eq(i)*alpha*q_eq + A_eq*eps*q_BL(3)) / (eps*sigma*A_eq)).^(1/4) - 273;
        Q_solar_rad_eq(i) = sin(theta(i))*q_eq;
    else
        Q_inf_e_op(i) = A_eq * (eps*sigma*T_op^4 - eps*q_BL(3)) - 20;
        Q_solar_rad_eq(i) = 0;
        T_e_op(i) = ((20 + A_eq*eps*q_BL(3)) / (eps*sigma*A_eq)).^(1/4) - 273;
    end
end

figure(3)
grid on
yyaxis left
plot(periodE,Q_inf_e_op,'LineWidth',3)
xlabel('Time (Hr)')
ylabel('Power (W)')
yyaxis right
plot(periodE,T_e_op,'LineWidth',3)
xlabel('Time (Hr)')
ylabel('Temperature (^{\circ}C)')
title('Operational Equinox with Eclipse')

%% survival heater power as a function of time

% designed to be held at -40 deg C

% winter solstice
for i=1:length(theta)
    if (theta(i) <= pi) % in sunlight
        Q_inf_w_su(i) = A_eq * (eps*sigma*T_su^4  - eps*q_BL(1)) - A_rot_eq(i)*alpha*q_w;
        T_w_su(i) = ((A_rot_eq(i)*alpha*q_w + A_eq*eps*q_BL(1)) / (eps*sigma*A_eq)).^(1/4) - 273;
    elseif (theta(i) > pi)
        Q_inf_w_su(i) = A_eq * (eps*sigma*T_su^4 - eps*q_BL(1));
        T_w_su(i) = ((A_eq*eps*q_BL(1)) / (eps*sigma*A_eq)).^(1/4) - 273;
    end
end

figure(4)
grid on
yyaxis left
plot(periodE,Q_inf_w_su,'LineWidth',3)
xlabel('Time (Hr)')
ylabel('Power (W)')
yyaxis right
plot(periodE,T_w_su,'LineWidth',3)
xlabel('Time (Hr)')
ylabel('Temperature (^{\circ}C)')
title('Survival Winter Solstice')

% summer solstice
for i=1:length(theta)
    if (theta(i) <= pi) % in sunlight
        Q_inf_s_su(i) = A_eq * (eps*sigma*T_su^4  - eps*q_BL(2)) - A_rot_eq(i)*alpha*q_s;
        T_s_su(i) = ((A_rot_eq(i)*alpha*q_s + A_eq*eps*q_BL(2)) / (eps*sigma*A_eq)).^(1/4) - 273;
    elseif (theta(i) > pi)
        Q_inf_s_su(i) = A_eq * (eps*sigma*T_su^4 - eps*q_BL(2));
        T_s_su(i) = ((A_eq*eps*q_BL(2)) / (eps*sigma*A_eq)).^(1/4) - 273;
    end
end

figure(5)
grid on
yyaxis left
plot(periodE,Q_inf_s_su,'LineWidth',3)
xlabel('Time (Hr)')
ylabel('Power (W)')
yyaxis right
plot(periodE,T_s_su,'LineWidth',3)
xlabel('Time (Hr)')
ylabel('Temperature (^{\circ}C)')
title('Survival Summer Solstice')

% equinox
for i=1:length(theta)
    if (theta(i) > (pi - thetaE/2)) && (theta(i) < (pi + thetaE/2)) % eclipse
        Q_inf_e_su(i) = A_eq * (sigma*eps*T_su^4 - eps*q_BL(4));
        T_e_su(i) = ((A_eq*eps*q_BL(4)) / (sigma*eps*A_eq)).^(1/4) - 273;
    elseif (theta(i) <= pi) % in sunlight
        Q_inf_e_su(i) = A_eq * (eps*sigma*T_su^4  - eps*q_BL(3)) - A_rot_eq(i)*alpha*q_eq;
        T_e_su(i) = ((A_rot_eq(i)*alpha*q_eq + A_eq*eps*q_BL(3)) / (eps*sigma*A_eq)).^(1/4) - 273;
    else
        Q_inf_e_su(i) = A_eq * (eps*sigma*T_su^4 - eps*q_BL(3));
        T_e_su(i) = ((A_eq*eps*q_BL(3)) / (eps*sigma*A_eq)).^(1/4) - 273;
    end
end

figure(6)
grid on
yyaxis left
plot(periodE,Q_inf_e_su,'LineWidth',3)
xlabel('Time (Hr)')
ylabel('Power (W)')
yyaxis right
plot(periodE,T_e_su,'LineWidth',3)
xlabel('Time (Hr)')
ylabel('Temperature (^{\circ}C)')
title('Survival Equinox with Eclipse')

figure(7)
plot(periodE,Q_solar_rad_w,'LineWidth',2)
hold on
plot(periodE,Q_solar_rad_s,'LineWidth',2)
plot(periodE,Q_solar_rad_eq,'LineWidth',2)
legend('Winter Solstice','Summer Solstice','Equinox')
hold off
grid on
xlabel('Time (Hr)')
ylabel('Solar Flux (W/m^2)')
title('Solar radiation flux exposed to radiator for one day')


