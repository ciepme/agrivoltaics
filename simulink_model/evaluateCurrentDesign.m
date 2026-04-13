clear;
clc;
% for running the model once at a specific design vector (from agri var def) and getting clear results without having to open up simulink 
agrivoltaics_variable_definition;
x = agriVarStruct2Array(agriVar, agriParams);
results = agrivoltaic_wrapper(x, agriParams);

E = results(1);
P = results(2);
social_cost = results(3);
pv_revenue = results(4);
crop_revenue = results(5);
yearly_biomass = results(6);
total_panels = results(7);
social_value = -social_cost;
emissions_value = 190 * (E / 1000);

fprintf('current design vector from agrivoltaics_variable_definition.m:\n');
disp(x);

fprintf('model evaluation results:\n');
fprintf('  Social value      : $%.4f\n', social_value);
fprintf('  Profit (P)        : $%.4f\n', P);
fprintf('  Emissions value   : $%.4f\n', emissions_value);
fprintf('  PV revenue        : $%.4f\n', pv_revenue);
fprintf('  Crop revenue      : $%.4f\n', crop_revenue);
fprintf('  Yearly biomass    : %.4f kg/yr\n', yearly_biomass);
fprintf('  Total panels      : %.1f\n', total_panels);
fprintf('  CO2 displaced     : %.4f tons\n', E / 1000);
