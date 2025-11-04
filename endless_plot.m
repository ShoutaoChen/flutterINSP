clear all; clc;
load("all_data_Hong_Kong.mat");
plot_probability_cloud(flutter_safety_matrix, d_tot_range, d_start_range, global_params);

function plot_probability_cloud(flutter_safety_matrix, d_tot_range, d_start_range, global_params)
    % Enhanced plotting function: Plot flutter safety probability cloud map
    % Adaptive color range, red for small values, blue for large values
    
    figure('Position', [100, 100, 450, 400]);
    
    % Create custom color map - from red to blue
    n_colors = 256;
    custom_colormap = create_custom_colormap(n_colors);
    
    % Plot cloud map
    imagesc(d_start_range, d_tot_range, flutter_safety_matrix);
    
    % Set color map and range
    colormap(custom_colormap);
    
    % Adaptive color range (keeping 5% boundary margin)
    data_min = min(flutter_safety_matrix(:));
    data_max = max(flutter_safety_matrix(:));
    data_range = data_max - data_min;
    caxis([data_min - 0.05*data_range, 1]);
    
    % Add color bar and set label
    c = colorbar;
    c.Label.String = 'Overall flutter occurrence probability';
    c.Label.FontSize = 10;
    c.Label.FontName = 'Times';
    
    
    % Axis labels
    xlabel('Start day of year (d)', 'FontSize', 12, 'FontName', 'Times');
    ylabel('Erection duration (d)', 'FontSize', 12, 'FontName', 'Times');
    
    % Axis settings
    axis xy;
    axis tight;
    
    % Grid and border enhancement
    grid on;
    set(gca, 'GridAlpha', 0.3, 'LineWidth', 1);
    box on;
    
    % Add contour lines - 0.8 with red dashed line, 0.9 with blue dashed line
    hold on;
    
    % 0.8 contour - red dashed line, show only one label
    [C1, h1] = contour(d_start_range, d_tot_range, flutter_safety_matrix, ...
                       [0.8, 0.8], 'LineWidth', 1, 'Color', 'red', 'LineStyle', '--', ...
                       'ShowText', 'off'); % First turn off automatic labels
    
    % Manually add one 0.8 label
    if ~isempty(C1)
        % Find starting point of first contour segment
        idx = 1;
        while idx < size(C1, 2)
            level = C1(1, idx);
            n_points = C1(2, idx);
            if level == 0.8
                x_pos = C1(1, idx+1);
                y_pos = C1(2, idx+1);
                text(x_pos, 0.95*y_pos, '0.8', 'FontSize', 12, 'FontWeight', 'bold', ...
                     'Color', 'black', 'FontName', 'Times');
                break;
            end
            idx = idx + n_points + 1;
        end
    end
    
    % 0.9 contour - blue dashed line, show only one label
    [C2, h2] = contour(d_start_range, d_tot_range, flutter_safety_matrix, ...
                       [0.9, 0.9], 'LineWidth', 1, 'Color', 'blue', 'LineStyle', '--', ...
                       'ShowText', 'off'); % First turn off automatic labels
    
    % Manually add one 0.9 label
    if ~isempty(C2)
        % Find starting point of first contour segment
        idx = 1;
        while idx < size(C2, 2)
            level = C2(1, idx);
            n_points = C2(2, idx);
            if level == 0.9
                x_pos = C2(1, idx+1);
                y_pos = C2(2, idx+1);
                text(x_pos, y_pos+10, '0.9', 'FontSize', 12, 'FontWeight', 'bold', ...
                     'Color', 'black', 'FontName', 'Times');
                break;
            end
            idx = idx + n_points + 1;
        end
    end
    
    hold off;
    

    % Set x-axis range and ticks
    xlim([1, 365]);
    ylim([30, 360]);
    xticks(1:60:365);  % One tick every 60 days
    xticklabels({'1', '60', '120', '180', '240', '300', '360'});
    yticks(30:30:360);  % One tick every 60 days
    % Set graphic properties
    set(gcf, 'Color', 'white');
    set(gca, 'FontSize', 14, 'FontName', 'Times');
    text(0.95, 0.05, global_params.city_name, ...
     'Units', 'normalized', ...
     'FontSize', 14, ...
     'FontName', 'Times New Roman', ...
     'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'bottom', ...
     'Margin', 3);
end

function custom_cmap = create_custom_colormap(n)
    % Create custom color map: red -> orange -> yellow -> cyan -> blue
    % Corresponding to low to high values: red -> orange -> yellow -> cyan -> blue
    
    % Define key color points
    colors = [0.7  0.0  0.0;   % Dark red (low values)
              0.9  0.2  0.1;   % Red
              1.0  0.5  0.1;   % Orange-red
              1.0  0.8  0.2;   % Orange
              0.8  0.9  0.4;   % Yellow-green
              0.6  0.9  0.8;   % Cyan-blue
              0.4  0.8  1.0;   % Light blue
              0.2  0.6  1.0;   % Blue
              0.0  0.0  0.8];  % Dark blue (high values)
    
    % Create interpolation positions
    positions = linspace(0, 1, size(colors, 1));
    
    % Interpolate to generate complete color map
    custom_cmap = interp1(positions, colors, linspace(0, 1, n));
    
    % Ensure color values are within [0,1] range
    custom_cmap = max(0, min(1, custom_cmap));
end