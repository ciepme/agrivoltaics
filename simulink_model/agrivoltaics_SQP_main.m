
clear;
clc;
rng default;

agrivoltaics_variable_definition;

x0 = agriVarStruct2Array(agriVar, agriParams); %initial guess pulled from variable definition file
use_scaled_SQP = true; % adding switch for running with scaling or without
scale = [1e2, 1e2, 1e2, 1e3, 1e5, 1e5, 1e5]; % scaling values derived from the Hessian at last x_star
x0_scaled = x0 .* scale;
lb_scaled = lb .* scale;
ub_scaled = ub .* scale;

% running SQP objective function with scaled DVs
obj_scaled = @(x_scaled) agrivoltaic_social_cost_of_carbon_wrapper(x_scaled ./ scale, agriParams);


options = optimoptions('fmincon','Algorithm', 'sqp','Display', 'iter','StepTolerance', 1e-6,'OptimalityTolerance', 1e-6);

tic; % Start timer
if use_scaled_SQP
    [x_opt_scaled, fval, exitflag, output] = fmincon(obj_scaled, x0_scaled, [], [], [], [], lb_scaled, ub_scaled, [], options);
    x_opt = x_opt_scaled ./ scale;
else
    [x_opt, fval, exitflag, output] = fmincon(@(x) agrivoltaic_social_cost_of_carbon_wrapper(x, agriParams), x0, [], [], [], [], lb, ub, [], options);
end
time_taken = toc;

% Results
fprintf('Exit Flag: %d\n', exitflag);
fprintf('Time Taken: %.2f seconds\n', time_taken);
fprintf('Maximized Social Value: $%.2f\n', -fval); % Flip sign back to positive
disp(' ');
disp('Optimal Variables Found:');
fprintf('  Panel Height (z_p) : %.4f m\n', x_opt(1));
fprintf('  Panel Length (l_p) : %.4f m\n', x_opt(2));
fprintf('  Panel Width (w_p)  : %.4f m\n', x_opt(3));
fprintf('  Azimuth (phi)      : %.4f rad (%.4f deg)\n', x_opt(4), rad2deg(x_opt(4)));
fprintf('  Tilt (sigma)       : %.4f rad (%.4f deg)\n', x_opt(5), rad2deg(x_opt(5)));
fprintf('  Row Spacing (y_p)  : %.6f m\n', x_opt(6));
fprintf('  Panel Gap (x_p)    : %.6f m\n', x_opt(7));
