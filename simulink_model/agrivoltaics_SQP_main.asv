
clear;
clc;
rng default;

agrivoltaics_variable_definition;

x0 = agriVarStruct2Array(agriVar); %initial guess pulled from variable definition file



options = optimoptions('fmincon','Algorithm', 'sqp','Display', 'iter','StepTolerance', 1e-6,'OptimalityTolerance', 1e-6);

tic; % Start timer
[x_opt, fval, exitflag, output] = fmincon(@agrivoltaic_social_cost_of_carbon_wrapper, x0, [], [], [], [], lb, ub, [], options);
time_taken = toc;

% Results
fprintf('Exit Flag: %d\n', exitflag);
fprintf('Time Taken: %.2f seconds\n', time_taken);
fprintf('Maximized Social Value: $%.2f\n', -fval); % Flip sign back to positive
disp(' ');
disp('Optimal Variables Found:');
fprintf('  Panel Height (z_p) : %.2f m\n', x_opt(1));
fprintf('  Panel Length (l_p) : %.2f m\n', x_opt(2));
fprintf('  Panel Width (w_p)  : %.2f m\n', x_opt(3));
fprintf('  Azimuth (phi)      : %.2f rad (%.1f deg)\n', x_opt(4), rad2deg(x_opt(4)));
fprintf('  Tilt (sigma)       : %.2f rad (%.1f deg)\n', x_opt(5), rad2deg(x_opt(5)));
fprintf('  Row Spacing (y_p)  : %.2f m\n', x_opt(6));
fprintf('  Panel Gap (x_p)    : %.2f m\n', x_opt(7));