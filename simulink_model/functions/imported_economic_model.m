function totalCost = imported_economic_model(tracking,efficiency,moduleheight,modulewidth,modulelength,netpower,projectduration)

%Module Cost Parameters (Same between Fixed and Single)
moduleyieldloss = 0.01; % $lost/$directcost measure of upstream modules produced but scrapped
moduleweight = 12.56; %kg/m2 representative of UPV c-Si modules
celltomodule = 0.98; %losses associated with cells per module
moduledepreciation = 0.118; % $/yr per $ invested; 45% equipment and 55% facility
moduleprofitmsp = 0.25; %minimum profit for most us companies

%Intrinsic unit Cell cost elements of monofacial 400Wdc c-Si module assembled in the US w/
%cells imported form SE Asia; domestic backsheet and junction boxes; other
%components imported from China
cell = 87.72; % $cost per m2 intrinsic unit, cells from SE Asia
frame = 1.62; %cost per m intrinsic unit; extruded aluminum frame from China
frontglass = 4.78; %cost per m2 intrinsic unit 3.2 mm glass imported from various countries
backsheet = 2.56; %cost per intrinsic unit m2 Polymeric backsheet from domestic supplier
encapsulant = 2.43; %two layers of EVA sourced from China
junctionboxes = 3.74; %domestically produced with foreign diodes 70% enclosure, 30% diodes
othermaterial = 4.41; % ribbons, edge seals, pottants, adhesives from china
directlabor = 5.54; %direct labor on the manufacturing lines and labor for movement of materials 
facilitiesoverhead = 0.88; %includes utilities and maintenance
businessoverhead = 2.95 + (28630000/1500000);%fixed; salaried personal, warranty, and insurance coverage
depreciation = 56.87; %facility capex per annual kWdc of nameplate production capacity
shipping = 0.32; %domestic shipping, including packaging factor

%Intrinsic unit Inverter cost elements for a 350 Wac microinverters using
%domestic components
pcbassemblies = 186; %printed circuit board with compoments mounted
enclosure = 48.45; %polymeric enclosure and any parts not mounted to PCB
directlaborinvert = 33.93; %direct labor
facilitiesoverheadinvert = 1.7; %includes utilities and maintenance
businessoverheadinvert = 2.7+(58200000/2000000);
depreciationinvert = 47.98; %capex per annual kWac of nameplate production capacity, split between equipment
shippinginvert = 1.09;%domestic shipping, including packaging for shipment

%Inverter Cost Parameters (Same between Fixed and Single)
inverterweight =0.935; %kginverter/kWacinverter
inverterdepreciation = 0.118; %same as module depreciation
inverteryieldloss = 0.01; %same as module yield loss
inverterprofitmsp = 0.25;%$/yr profit/ $invested

%Other (Same between Fixed and Single)
salestaxrate = 0.058; %$tax/$cost
contingencyrate = 0.0215; % $reserved/$cost money reserved in case of delays/ obstacles in implementation
overheadrate = 0.01; % $ overhead/ $ cost cost of business management 

%Structural Balance of System (SBOS) Parameters Single Axis
SingleAxis_SBOS_FastenerWeight = 0.2; %kg/m2 modules
SingleAxis_SBOS_TorqueTubeWeight = 4.8; %kg/m2 modules
SingleAxis_SBOS_EBOSPadCapacity = 50; %kWac/m2 pad

%Electrical Balance of System (EBOS) Parmeter Single Axis
SingleAxis_EBOS_EBOSweight = 2; %kgEBOS/kWac inverter

%System Economic Parameters Single Axis
SingleAxis_ILR = 1.34; %kWdc modules/ kW ac inverter: Inverter loading ratio for tracking UPV systems
SingleAxis_Tariff301 = 0.25; %$ tariff / $ cost, tariff costs applicable per 2025 Q1

%Operations and Maintenance Parameters Single Axis
SingleAxis_ModuleLossRate = 0.002; %fraction/year of modules lost
SingleAxis_InverterLossRate = 0.09; %fraction/year of inverters lost
SingleAxis_SBOSLossRate = 0.002; %fraction/year of SBOS lost
SingleAxis_EBOSLossRate = 0.002; %fraction/year of EBOS lost
SingleAxis_InsuranceRate = 0.0025; %$premium/year / $cost

%%Calculations

%ModuleCostCalculation
ModuleCostWithoutProfit = cell*((1/celltomodule)*(1/(1-moduleyieldloss)))+frame*(2*(modulewidth+modulelength)/(modulewidth*modulelength*efficiency)*(1/(1-moduleyieldloss)))+ ...
    (frontglass+backsheet+encapsulant)*((1/efficiency)*(1/(1-moduleyieldloss)))+(junctionboxes+othermaterial)*(1/(efficiency*modulelength*modulewidth)*(1/(1-moduleyieldloss)))+...
    (directlabor*facilitiesoverhead)*(1/efficiency)+businessoverhead+depreciation*moduledepreciation+shipping*(moduleweight/efficiency);
ModuleCost = moduleprofitmsp*(ModuleCostWithoutProfit-depreciation*moduledepreciation+depreciation)/4+ModuleCostWithoutProfit;

%InverterCostCalculation
InverterCostWithoutProfit = (enclosure+pcbassemblies)*(1/(1-inverteryieldloss))+directlaborinvert+facilitiesoverheadinvert+businessoverheadinvert;
InverterCost = inverterprofitmsp*(InverterCostWithoutProfit+depreciation)/4+InverterCostWithoutProfit+depreciationinvert*inverterdepreciation+shippinginvert*inverterweight;

if tracking == 0
    heightdelta = moduleheight-2.5; %set the number subtracted here to be the lower bound of the system; this bound was chosen to allow for minimum tractor clearance
    rackingheightelasticity = 0.35; %height driven scaling factor; adjustable, based on NREL data
    freightweightelasticity = 0.23; %freight adjusted scaling factor for increased weight due to height of module
    Piersweight = 5.68*(1+rackingheightelasticity*heightdelta^2);
    Shippingcosting = 7*(1+freightweightelasticity*heightdelta);
    SBOS = Piersweight+2.32*1.5+1.07+0.36+346.4*efficiency/(SingleAxis_ILR*SingleAxis_SBOS_EBOSPadCapacity)+Shippingcosting;
    EBOS = (30.01+5.46+6.71)*SingleAxis_ILR+51.08+4.69+12.03+56.59+51.26+2.26+(15.14+1883000/74627)+131000/74627;
    laborcomplexityfactor = 1+(0.25*heightdelta);
    rentaldurationfactor = 1+(0.3*heightdelta);
    Fieldwork = (9.48*laborcomplexityfactor)+52.86*1.2+(1+0.1*heightdelta)*(10.54)+6.84+(9*rentaldurationfactor)+1.39;
    Officework=(SingleAxis_ILR/efficiency)*(1.04+0.1+1.82)+(0.21+58.39+3.31)+(131200+478600+411400+32800)/74627;
    Salestax = (ModuleCost+InverterCost+SBOS+EBOS)*salestaxrate;
    Contingency = (ModuleCost+InverterCost+SBOS+EBOS+Fieldwork+Officework)*contingencyrate;
    Management = (ModuleCost+InverterCost+SBOS+EBOS+Fieldwork+Officework)*overheadrate+1363000/200000;
    Other=Salestax+Contingency+Management;
    OM = (0.55+0.65)*(1/efficiency)+ModuleCost*SingleAxis_ModuleLossRate+InverterCost*SingleAxis_InverterLossRate/SingleAxis_ILR+SBOS*SingleAxis_SBOSLossRate/efficiency+...
        EBOS*SingleAxis_EBOSLossRate/SingleAxis_ILR+SingleAxis_InsuranceRate*(SBOS+EBOS+ModuleCost+InverterCost);
else
    heightdelta = moduleheight-2.5; %set the number subtracted here to be the lower bound of the system; this bound was chosen to allow for minimum tractor clearance
    rackingheightelasticity = 0.35; %height driven scaling factor; adjustable, based on NREL data
    freightweightelasticity = 0.2; %freight adjusted scaling factor for increased weight due to height of module
    Piersweight = 5.68*(1+rackingheightelasticity*heightdelta^2);
    Shippingcosting = 10*(1+freightweightelasticity*heightdelta);
    SBOS = 7.44+Piersweight+2.32+1.07+0.36+3.63+3.28+1.83+1.42+346.4*efficiency/(SingleAxis_ILR*SingleAxis_SBOS_EBOSPadCapacity)+Shippingcosting;
    EBOS = (30.01+5.46+6.71)*SingleAxis_ILR+51.08+4.69+12.03+56.59+51.26+2.26+(15.14+1883000/74627)+131000/74627;
    laborcomplexityfactor = 1+(0.25*heightdelta);
    rentaldurationfactor = 1+(0.3*heightdelta);
    Fieldwork = (9.48*laborcomplexityfactor)+52.86*(efficiency/SingleAxis_ILR)+(1+0.1*heightdelta)*(6.92+1654000/456621)+6.84+(11.48*rentaldurationfactor)+1.39;
    Officework=(SingleAxis_ILR/efficiency)*(1.04+0.1+1.82)+(0.21+58.39+3.31)+(131200+478600+411400+32800)/74627;
    Salestax = (ModuleCost+InverterCost+SBOS+EBOS)*salestaxrate;
    Contingency = (ModuleCost+InverterCost+SBOS+EBOS+Fieldwork+Officework)*contingencyrate;
    Management = (ModuleCost+InverterCost+SBOS+EBOS+Fieldwork+Officework)*overheadrate+1363000/200000;
    Other=Salestax+Contingency+Management;
    OM = (0.55+0.65)*(1/efficiency)+ModuleCost*SingleAxis_ModuleLossRate+InverterCost*SingleAxis_InverterLossRate/SingleAxis_ILR+SBOS*SingleAxis_SBOSLossRate/efficiency+...
        EBOS*SingleAxis_EBOSLossRate/SingleAxis_ILR+SingleAxis_InsuranceRate*(SBOS+EBOS+ModuleCost+InverterCost);
end
CostMSP = ModuleCost+(InverterCost+EBOS+Officework)*(1/SingleAxis_ILR)+(SBOS+Fieldwork)*(1/efficiency)+Other+OM*(projectduration/20);
totalCost = netpower*CostMSP;
