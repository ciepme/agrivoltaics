function P = pv_model(DNI, DHI, beta_s, phi_s, l_p, w_p, phi, sigma, n_p)

%PV model calculating incidence angle, total radiation, and estimated power
%output
arguments (Input)
    DNI ; %direct normal irradiance
    DHI ; % diffuse horizontal irradiance
    beta_s ; %sun altitude angle (angle of sun above the local horizon)
    phi_s ; %sun azimuth angle (angle between the sun and true south)
    l_p ; % panel length (m)
    w_p ; % panel width (m)
    phi ; % azimuth angle (radians) - relative to true South, going ccw e.g. pi/2 rad is East
    sigma ; % tilt angle (radians) - fixed sloping angle of PV relative to horizontal (xy) plane 
    n_p; %panel efficiency
end

arguments (Output)
    P ;%estimated power output per panel
end

p_tot=0;
for i = 1:24
    if beta_s(i) <= 0
        continue; % Skip to next hour (P_inst remains 0 effectively)
    end
% Incidence angle
cos_theta = cos(beta_s(i))*cos(phi_s(i) - phi)*sin(sigma) ...
          + sin(beta_s(i))*cos(sigma);

theta = acos(cos_theta);

% Direct beam
I_db = DNI(i) .* cos_theta;

% Diffuse irradiance
I_diff = DHI(i) .* (1 + cos(sigma)) / 2;

% Reflected (Need to double check equation)
%I_ref = (DNI*sin(beta_s) + DHI) .* albedo .* (1 - cos(tilt)) / 2;

% Total Intercepted Radiation
R = I_db + I_diff; %+ I_ref;

% Simple power model
A_p = l_p*w_p ; %panel area m^2
p_hour = n_p * A_p * R; % power equation
p_tot=p_tot+p_hour;
end
P=p_tot;
end