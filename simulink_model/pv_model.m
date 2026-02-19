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

% Incidence angle
cos_theta = cos(beta_s(1))*cos(phi_s(1) - phi)*sin(sigma) ...
          + sin(beta_s(1))*cos(sigma);

theta = acos(cos_theta);

% Direct beam
I_db = DNI(1) .* cos_theta;

% Diffuse irradiance
I_diff = DHI(1) .* (1 + cos(sigma)) / 2;

% Reflected (Need to double check equation)
%I_ref = (DNI*sin(beta_s) + DHI) .* albedo .* (1 - cos(tilt)) / 2;

% Total Intercepted Radiation
R = I_db + I_diff; %+ I_ref;

% Simple power model
A_p = l_p*w_p ; %panel area m^2
P = n_p * A_p * R; % power equation

end
