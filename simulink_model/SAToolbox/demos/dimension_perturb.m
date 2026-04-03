function [d_new, t_new] = dimension_perturb(d_old,t_old)

d_local_mod = 0.1;
t_local_mod = 0.05;

d_new = d_old + d_local_mod .* (2.*(rand(1) - 0.5));
t_new = t_old + t_local_mod .* (2.*(rand(1) - 0.5));