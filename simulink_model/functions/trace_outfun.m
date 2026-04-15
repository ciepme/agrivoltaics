function stop = trace_outfun(~, optimValues, state)
% Records per iteration convergence data into base workspace as
% conv_history.mat--> this function is called after every SQP  iteration to
% record the data below:
% columns: [iteration,social_value, first_order_optimality, step_size]
% not sure what to do with first_order_optimality and step_size values yet, but
% might as well record them
    persistent data;
    if strcmp(state, 'init')
        data = [];
    elseif strcmp(state, 'iter')
        data(end+1,:) = [optimValues.iteration, ...
                         -optimValues.fval, ...       % social value-->positive
                         optimValues.firstorderopt, ... 
                         optimValues.stepsize];
        assignin('base', 'conv_history', data);
    end
    stop = false;
end