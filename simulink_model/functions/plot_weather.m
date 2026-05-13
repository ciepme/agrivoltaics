% Load the seasonal weather data we saved previously
load('pv_weather_4seasons.mat', 'weather_struct');

% Define the seasons to loop through
season_names = {'spring', 'summer', 'fall', 'winter'};
plot_titles = {'Spring Average', 'Summer Average', 'Fall Average', 'Winter Average'};

% Create a new figure
figure('Name', 'Seasonal Average DNI and DHI', 'Position', [100, 100, 900, 600]);

for i = 1:4
    curr_season = season_names{i};
    
    % Extract the DNI and DHI arrays
    dni = weather_struct.(curr_season).DNI;
    dhi = weather_struct.(curr_season).DHI;
    
    % Create an x-axis for the hours (assuming 24 hourly data points)
    hours = 0:(length(dni)-1); 
    
    % Create a 2x2 subplot
    subplot(2, 2, i);
    
    % Plot DNI and DHI
    plot(hours, dni, '-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'DNI (Direct)');
    hold on;
    plot(hours, dhi, '-s', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'DHI (Diffuse)');
    hold off;
    
    % Formatting
    title(plot_titles{i}, 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Hour of Day');
    ylabel('Irradiance (W/m^2)');
    xlim([0 23]);
    
    % Keep the y-axis consistent across all plots for easy visual comparison
    ylim([0 max(max(weather_struct.summer.DNI), 1000) * 1.1]); 
    
    grid on;
    legend('Location', 'northeast');
end

sgtitle('Average Daily Hourly Irradiance by Season (Ann Arbor, MI)', 'FontSize', 14, 'FontWeight', 'bold');