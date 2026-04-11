clear;
clc;

agrivoltaics_variable_definition;

% fixed axis SQP design vector:
%           [z_p, l_p, w_p, phi, sigma, y_p, x_p]
x_star = [4.5, 2.5, 1.5, deg2rad(17.5),deg2rad(37.0), 2.5488, 0.1129];
% x_star = [3.2, 2.0, 1.1, deg2rad(17.5),deg2rad(37.0), 2.5488, 0.1129];
% x_star = [4.4, 2.4, 1.4, deg2rad(17.5),deg2rad(37.0), 2.5488, 0.1129];

names = {'PV_z_p', 'PV_l_p', 'PV_w_p', 'PV_phi', 'PV_sigma', 'PV_y_p', 'PV_x_p'};

% calling wrapper to minimize social cost/ maximize social value
obj = @(x) agrivoltaic_social_cost_of_carbon_wrapper(x, agriParams);

range = ub - lb; % bounds of each DV
h = 1e-3 .* range; % % finite difference step size is based off of the constrained range of each DV
f0 = obj(x_star); % model evalues at x star

n = numel(x_star); % 7 DVs
Hdiag = zeros(n,1); % for storing the Hessian diagonal entries
% finite difference uses different method depending on if constraint of DV
% is tight (backwards/normal/forwards)
finiteDiffMethod = cell(n,1);


for i = 1:n
    ei = zeros(1,n); 
    ei(i) = 1;

    % using central diff method
    if x_star(i) - h(i) >= lb(i) && x_star(i) + h(i) <= ub(i)
        fp = obj(x_star + h(i) .* ei); %increasing var i val by h(i), then evaluating at that design point
        fm = obj(x_star - h(i) .* ei);
        Hdiag(i) = (fp - 2 .* f0 + fm) ./ (h(i).^2); % central diff formula to approx. second deriv wrt var i
        finiteDiffMethod{i} = 'central';
    
    % if DV is at upper constraint
    elseif x_star(i) - 2 .* h(i) >= lb(i) %checking if DV can move backward 2 steps w/out violating constraint
        f1 = obj(x_star - h(i) .* ei); % evaluating one step back
        f2 = obj(x_star - 2 .* h(i) .* ei);
        f3 = obj(x_star - 3 .* h(i) .* ei); % 3 steps of h(i) backward
        Hdiag(i) = (f0 - 2 .* f1 + f2) ./ (h(i).^2); % standard second order backward finite diff formula
        finiteDiffMethod{i} = 'backward';

    elseif x_star(i) + 2 .* h(i) <= ub(i)
        f1 = obj(x_star + h(i) .* ei);
        f2 = obj(x_star + 2 .* h(i) .* ei);
        f3 = obj(x_star + 3 .* h(i) .* ei);
        Hdiag(i) = (f0 - 2 .* f1 + f2) ./ (h(i).^2); % forward diff formula
        finiteDiffMethod{i} = 'forward';

    else
        error('Finite Diff method for Hessian did not work')
    end
end

results = table(names', x_star', lb', ub', h', finiteDiffMethod, Hdiag, ...
    'VariableNames', {'variable','x_star','lb','ub','step','finiteDiffMethod','H_ii'});

disp(results);
fprintf('\nThese H_ii values are for the minimized objective (social_cost).\n');
fprintf('For the Hessian of social value, multiply all entries by -1.\n');
