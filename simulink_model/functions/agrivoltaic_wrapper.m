function results = agrivoltaic_wrapper(custom_var, agriParams)

    %set custom variables
    agriVar = agriVarArray2Struct(custom_var, agriParams);
    agriVar = orderfields(agriVar);

    %create Simulation objects
    agriSim = Simulink.SimulationInput('agrivoltaics_v1');

    %set local as default variable (normal looks in base workspace)
    agriSim = agriSim.setVariable('agriVar', agriVar);
    agriSim = agriSim.setVariable('agriParams', agriParams);

    %set in rapid accelerator
    %agriSim = agriSim.setModelParameter('SimulationMode','Rapid');
    %agriSim = agriSim.setModelParameter('RapidAcceleratorUpToDateCheck','off'); 

    ensure_worker_bus_objects(agriParams, agriVar);

    %run sim
    out = sim(agriSim);

    %fprintf('P: %.2f\n', out.p.Data);

    E = out.e.Data; % total CO2 displaced via solar power generation across 1 year
    P = out.p.Data; % profit

    social_cost = -1.*(P + 190 .* (E ./ 1000));

    crop_revenue = out.crop_revenue.Data;
    total_panels = out.total_panels.Data;
    pv_revenue = out.pv_revenue.Data;
    yearly_biomass = out.yearly_biomass.Data;
    yearly_energy = out.yearly_energy.Data;

    results = [E, P, social_cost, pv_revenue, crop_revenue, yearly_biomass, total_panels, yearly_energy];
    closeExcel;
end

function ensure_worker_bus_objects(agriParams, agriVar)
    signature = bus_signature(agriParams, agriVar);

    has_params_bus = evalin('base', 'exist(''params_bus'', ''var'') == 1');
    has_var_bus = evalin('base', 'exist(''var_bus'', ''var'') == 1');
    has_signature = evalin('base', 'exist(''agri_bus_signature'', ''var'') == 1');

    if has_params_bus && has_var_bus && has_signature
        existing_signature = evalin('base', 'agri_bus_signature');
        if isequal(existing_signature, signature)
            return;
        end
    end

    params_info = Simulink.Bus.createObject(orderfields(agriParams));
    params_bus = evalin('base', params_info.busName);
    assignin('base', 'params_bus', params_bus);
    evalin('base', sprintf('clear(''%s'')', params_info.busName));

    var_info = Simulink.Bus.createObject(orderfields(agriVar));
    var_bus = evalin('base', var_info.busName);
    assignin('base', 'var_bus', var_bus);
    evalin('base', sprintf('clear(''%s'')', var_info.busName));

    assignin('base', 'agri_bus_signature', signature);
end

function signature = bus_signature(agriParams, agriVar)
    signature = struct();
    signature.params_fields = field_signature(orderfields(agriParams));
    signature.var_fields = field_signature(orderfields(agriVar));
end

function sig = field_signature(value)
    names = fieldnames(value);
    sig = cell(numel(names), 4);

    for i = 1:numel(names)
        field_value = value.(names{i});
        sig{i, 1} = names{i};
        sig{i, 2} = class(field_value);
        sig{i, 3} = size(field_value);
        if isstruct(field_value)
            sig{i, 4} = field_signature(orderfields(field_value));
        else
            sig{i, 4} = [];
        end
    end
end
