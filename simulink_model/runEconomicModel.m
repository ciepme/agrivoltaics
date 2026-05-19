function Cost = runEconomicModel(tracking, efficiency, height, width, length, power)
% choosing the right txt files (either fixed or single)

override = 1
% PC Function (Excel economics)    
if 
    if (tracking == 0)
            targettxt = 'FixedAxis.txt'; 
        else
            targettxt = 'SingleAxis.txt';
        end
    % Create the absolute path to the file inside the subfolder
    %fileName = fullfile(subfolderPath, targettxt);
    fileName = fullfile(pwd, 'functions', 'PVeconomicsfiles', targettxt);
    excelFile = fullfile(pwd, 'functions', 'PVeconomicsfiles', 'EconomicModel.xlsm');
    persistent app wb excelFileCached
    
    if isempty(app) || ~iscom(app)
        try
            % Try to grab an existing Excel instance first
            app = actxGetRunningServer('Excel.Application');
        catch
            % If none, start a new one
            app = actxserver('Excel.Application');
        end
        app.DisplayAlerts = false;
        app.Visible = false; % Keep it hidden for speed
    end
    
    % Safety Check
    if ~exist(fileName, 'file')
        error('File not found: %s. Please check if the PVeconomicsfiles folder exists!', fileName);
    end
    
        try
            % Adjusting variables in the txt file
            content = fileread(fileName);
            newWidth = mean([width, length]);
            
            % Formatting the lines
            sizeinput = sprintf('System.Size=%.2f', power);
            effIn = sprintf('Factors.ModuleEfficiency=Module\\\\m2 modules\\\\kWdc modules\\\\%.3f\\\\', efficiency);
            hIn = sprintf('Factors.ModuleHeight=Module\\\\\\\\m\\\\%.3f\\\\\\\\', height);
            wIn = sprintf('Factors.ModuleWidth=Module\\\\\\\\m\\\\%.3f\\\\\\\\', newWidth);
    
            % Overwriting the original txt file
            updated = regexprep(content, '(?m)^System\.Size=[\d\.]+', sizeinput);
            updated = regexprep(updated, '(?m)^Factors\.ModuleEfficiency=[^\r\n]*', effIn);
            updated = regexprep(updated, '(?m)^Factors\.ModuleHeight=[^\r\n]*', hIn);
            updated = regexprep(updated, '(?m)^Factors\.ModuleWidth=[^\r\n]*', wIn);
    
            % saving new txt file over original
            fid = fopen(fileName, 'w');
            fprintf(fid, '%s', updated);
            fclose(fid);
    
            % pulling and inputing new txt file into excel model
            
            app.DisplayAlerts = false;
            app.Visible=false;
        if isempty(wb) || ~iscom(wb)
            wb = app.Workbooks.Open(excelFile);
        end
        % Robust Macro Name generation
    [~, wbName, wbExt] = fileparts(excelFile);
    
    % The format: "'WorkbookName.xlsm'!MacroName"
    % Note the single quotes inside the double quotes
    targetMacro = sprintf("'%s%s'!FileHandling.LoadFile", wbName, wbExt);
    
    % EXECUTE: Pass the 'fileName' (the .txt path) as the second argument
    app.Run(targetMacro, fileName);
        
        app.CalculateFull;
        while ~strcmp(app.CalculationState, 'xlDone')
            pause(0.1); % Wait for background calculations to finish
        end
        Cost = wb.Sheets.Item('System').Range('H22').Value;
        catch ME
            fprintf('Error: %s\n', ME.message);
            Cost = NaN;
        end
else 

        fixedaxis_agrivoltaic_capital_cost_rate = 1500;
        singleaxis_agrivoltaic_capital_cost_rate = 2500;
        if tracking == 0 
            Cost=power*fixedaxis_agrivoltaic_capital_cost_rate;
        else
            Cost=power*singleaxis_agrivoltaic_capital_cost_rate;
        end
end