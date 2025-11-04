clear all;clc;
%% Input module
cal_series = [1,2,3,4,5,6]; % Input city numbers to calculate simultaneously
                            % 1:Dalian, 2:Qingdao, 3:Hangzhou, 4:Taipei, 5:Xiamen, 6:Hong Kong
%Bridge to be analysed
bridge_type = 1; % 1-the Severn Bridge; 2-the Xihoumen Bridge
% Construction period range for calculation
d_tot_range = 30:10:360;    % Total erection duration
d_start_range = 1:1:365;   % Erection start date
%% Main function
folder_path = 'WDSP_database';
switch bridge_type
    case 1
        span = 988; % main span of bridge (m)
        z = 40;     % reference height at bridge deck (m)
        global_params.U_f = [59.04,46.29,47.94,55.39,61.22,64.95,65.81,64.59,62.21,59.22]; % Flutter critical wind speeds for each stage (Severn)
    case 2
        span = 1650;
        z = 50;
        global_params.U_f = [70.00,68.43,70.11,66.05,63.24,59.31,60.72,60.36,58.93,90.44]; % Flutter critical wind speeds for each stage (Xihoumen)
end
for i_seri = 1:length(cal_series)
    city_number = cal_series(i_seri); 
    city_names = {'Dalian', 'Qingdao', 'Hangzhou', 'Taipei', 'Xiamen', 'Hong Kong'};
    file_names = {'wind_dalian_30.txt', 'wind_qingdao_30.txt', 'wind_hangzhou_30.txt', 'wind_taibei_30.txt', 'wind_xiamen_30.txt', 'wind_hongkong_30.txt'};
    
    if city_number >= 1 && city_number <= 6
        wind_file = file_names{city_number};
        global_params.city_name = city_names{city_number};
    else
        error('Invalid city number, please enter a number between 1-6');
    end
    % Calculation of K, refer to "Wind-Resistant Design Specification for Highway Bridges" 4.2.6. Selected cities belong to risk area R1;
    surface = 'A';% Bridge site terrain roughness category A,B,C,D
    k_f=1.05;
    k_t =1;% Terrain condition coefficient
    k_h = cal_k_h(z, surface);
    % Calculation of Γ
    gamma_t = interpolate_gamma_t(span, surface, 'GAMMA_T.CSV');
    gamma_f = 1.15;  % Flutter stability partial factor, taken as 1.4 when using specification calculation for flutter critical wind speed,
                 % taken as 1.15 when using wind tunnel test method to obtain flutter critical wind speed,
                % taken as 1.25 when using virtual wind tunnel test method;
    gamma_a = 1.0;
   
    global_params.K = k_t*k_f*k_h;        % Coefficient K
    global_params.Gamma = gamma_t*gamma_f*gamma_a;    % Coefficient Gamma  
    global_params.N = length(global_params.U_f);          % Number of construction stages
    
    % Read data
    data = read_wind_data(wind_file,folder_path);
    % Terrain correction
    U_corrected = terrain_correction(data,city_number);
    
    % Generate feasible region cloud map
    % Input:
    %   global_params - Global parameters

    
    % Initialize result matrices
    flutter_safety_matrix = zeros(length(d_tot_range), length(d_start_range));
    P_so_matrix = zeros(length(d_tot_range), length(d_start_range));
    
    % Progress display
    h = waitbar(0, 'Calculating...');
    total_iterations = length(d_tot_range) * length(d_start_range);
    current_iteration = 0;
    
    % Traverse all combinations
    for i = 1:length(d_tot_range)
        d_tot = d_tot_range(i);
        
        for j = 1:length(d_start_range)
            d_start = d_start_range(j);
            
            % Check flutter safety probability of single scheme
            [P_so_total, P_sf_stages] = check_single_scheme_probability(data,U_corrected,d_tot, d_start, global_params);
            
            flutter_safety_matrix(i, j) = P_so_total;
            P_so_matrix(i, j) = P_so_total;
            
            % Update progress bar
            current_iteration = current_iteration + 1;
            if mod(current_iteration, 1000) == 0
                waitbar(current_iteration/total_iterations, h);
            end
        end
    end
    close(h);
    
    % Plot flutter safety probability cloud map
    plot_probability_cloud(flutter_safety_matrix, d_tot_range, d_start_range, global_params);
    
    % === New: Save all data to MAT file ===
    city_name_for_file = strrep(global_params.city_name, ' ', '_'); % Replace spaces with underscores
    filename = sprintf('S_all_data_%s.mat', city_name_for_file);
    
    % Save all relevant variables
    save(filename, ...
        'city_number', 'global_params', 'data', 'U_corrected', ...
        'd_tot_range', 'd_start_range', 'flutter_safety_matrix', ...
        'P_so_matrix', 'span', 'surface', 'z', 'k_f', 'k_t', 'k_h', ...
        'gamma_t', 'gamma_f', 'gamma_a');
    
    fprintf('All data saved to file: %s\n', filename);
    % === End of new part ===
end

%% Sub functions
function [P_so_total, P_fs_stages] = check_single_scheme_probability(data,U_corrected,d_tot, d_start, global_params)
    % Modified function: Calculate flutter safety probability of single scheme
    % Input:
    %   d_tot - Total erection duration
    %   d_start - Erection start date
    %   global_params - Global parameters
    % Output:
    %   P_so_total - Overall flutter safety probability
    %   P_fs_stages - Flutter occurrence probability for each stage
    
    N = global_params.N;
    U_f = global_params.U_f;
    K = global_params.K;
    Gamma = global_params.Gamma;
    U_sf = U_f/K/Gamma;
    % Duration of each stage
    d = d_tot / N;
    
    P_fs_stages = zeros(1, N); % Flutter occurrence probability for each stage
    
    % Calculate weight coefficients w(i)
    w = (1 ./ U_f) / mean(1 ./ U_f);
    
    % Check each stage
    for i = 1:N
        % Calculate stage start time
        d_s_i = d_start + (i-1) * d;
        % Normalize to 1-365
        d_s_i = mod(d_s_i - 1, 365) + 1;

        % Calculate Gumbel distribution parameters for this stage
        [mu, alpha] = GM_params(data, U_corrected, d_s_i, d);
      
        % Calculate flutter occurrence probability P_sf for this stage
        P_fs_stages(i) = exp(-exp(-(U_sf(i)-mu)/(1.3*alpha)));
    end
    
    % Calculate overall flutter safety probability P_so = ∏(1 - P_sf(i))^w(i)
    % Or according to your description: P_so is the product of P_sf(i)^w(i) for each stage
    P_so_total = prod((P_fs_stages) .^ w);
end


function [mu, alpha] = GM_params(data, U_corrected, d_str, d_rp)
% Calculate Gumbel distribution parameters for the first stage
    EWS_history = group_first_period(data, U_corrected, d_str, d_rp);
    % Get all annual maximum wind speeds for the first stage (excluding 0 and NaN values)
    valid_data = EWS_history(EWS_history > 0 & ~isnan(EWS_history));
    
    if length(valid_data) < 2
        error('Insufficient valid data in first stage, unable to calculate reliable parameters');
    end
    
    % Calculate mean and standard deviation
    mean_EWS = mean(valid_data);
    std_EWS = std(valid_data);
    
    % Calculate Gumbel parameters
    alpha = sqrt(6) / pi * std_EWS;
    mu = mean_EWS - 0.5772 * alpha;
end

function data = read_wind_data(filename,floder_path)
% Read wind speed data file (unchanged)
    full_filename = fullfile(folder_path, filename);
    
    fprintf('Reading data file: %s\n', full_filename);
    
    % Read text file
    raw_data = load(filename);
    
    % Extract each column data
    data.station = raw_data(:, 1);    % Station number
    data.year = raw_data(:, 2);       % Year
    data.month = raw_data(:, 3);      % Month
    data.day = raw_data(:, 4);        % Day
    data.V_600 = raw_data(:, 5)*0.5144;      % V_600 (10min daily maximum wind speed,unit: m/s) original unit is knot
    data.V_3 = raw_data(:, 6)*0.5144;        % V_3 (3s daily maximum wind speed,unit: m/s) original unit is knot
    
    fprintf('Data reading completed, total %d records\n', length(data.year));
end

function U_corrected = terrain_correction(data,city_number)
% Terrain correction calculation (unchanged)
    % Initialize output as original data
    U_corrected = data.V_600;
    
    % Get all years
    years = unique(data.year);
    n_years = length(years);
    
    % Store annual average G0 values
    yearly_G0 = zeros(1, n_years);
    valid_years = false(1, n_years);
    
    % Calculate annual average G0 values
    for i = 1:n_years
        current_year = years(i);
        year_idx = (data.year == current_year);
        
        % Filter valid data: V_3 and V_600 within reasonable range
        valid_idx = year_idx & (data.V_3 > 0) & (data.V_3 < 100) & ...
                   (data.V_600 > 5) & (data.V_600 < 50);
        
        % Calculate G values
        G_values = data.V_3(valid_idx) ./ data.V_600(valid_idx);
        
        % Check valid data days (excluding abnormal G values)
        valid_G = G_values(G_values > 0.1 & G_values < 2);
        
        if length(valid_G) >= 3  % At least 3 days of valid data per year required
            yearly_G0(i) = mean(valid_G);
            valid_years(i) = true;
        else
            yearly_G0(i) = NaN;
            valid_years(i) = false;
        end
    end
    
    % Check if sufficient years for linear fitting
    valid_indices = find(valid_years);
    if length(valid_indices) < 3
        fprintf('Insufficient valid years (%d), terrain correction not performed\n', 3);
        U_corrected(data.V_600 >= 50 | data.V_600 <= 0) = NaN;
        return;
    end
    
    % Extract valid years and corresponding G0 values
    valid_years_vec = years(valid_years);
    valid_G0 = yearly_G0(valid_years);
    
    % Linear fitting G(t) = a*t + b
    t = valid_years_vec - min(valid_years_vec);  % Using minimum year as baseline

    p = polyfit(t, valid_G0, 1);

    % Calculate G(t) for all years
    all_t = years - min(valid_years_vec);
    G_t = polyval(p, all_t);
    

    % Perform terrain correction for each year
    for i = 1:n_years
        year_idx = (data.year == years(i));
        
        % Filter data requiring correction
        valid_correction = year_idx & (data.V_600 > 0) & (data.V_600 < 80);
        
        % Calculate z0 and β
        z0 = 10 / exp(2.32 / (G_t(i) - 1.08));
        if z0>0.5
            z0 = 0.5;
        end
        beta = (0.19 * (z0/0.05).^0.07 * log(10/z0))^(-1);
        
        % Apply correction
        U_corrected(valid_correction) = beta .* data.V_600(valid_correction);

    end
    
    % Process invalid data
    U_corrected(data.V_600 >= 50 | data.V_600 <= 0) = NaN;
    
    fprintf('Terrain correction completed, valid years %d/%d, fitting parameters: G(t)=%.4f*t+%.4f\n', ...
            length(valid_indices), n_years, p(1), p(2));
end

function EWS_history = group_first_period(data, U_corrected, d_str, d_rp)
% Calculate only annual maximum wind speed for the first stage (unchanged)
    years = unique(data.year);
    n_years = length(years);
    
    % Calculate total days in a year (based on actual data)
    first_year_data = data.year == years(1);
    days_in_year = 365;
    
    % Initialize output vector
    EWS_history = zeros(1, n_years);
    
    for i = 1:n_years
        year_data = data.year == years(i);
        year_U = U_corrected(year_data);
        year_day = calculate_day_of_year(data.month(year_data), data.day(year_data));

        % Calculate start and end days of first stage
        start_day = d_str;
        end_day = mod(d_str + d_rp - 1, days_in_year);
        
        % Determine day range included in stage
        if start_day <= end_day
            % Non-cross-year situation
            stage_idx = (year_day >= start_day) & (year_day <= end_day);
        else
            % Cross-year situation: from start_day to year-end, then from year-start to end_day
            stage_idx = (year_day >= start_day) | (year_day <= end_day);
        end
        
        % Extract wind speed data for this stage (excluding NaN values)
        stage_U = year_U(stage_idx);
        valid_stage_U = stage_U(~isnan(stage_U));
        
        if ~isempty(valid_stage_U)
            EWS_history(i) = max(valid_stage_U);
        else
            EWS_history(i) = NaN; % Set as NaN if no valid data
        end
    end
end

function day_of_year = calculate_day_of_year(month, day)
% Calculate day of year (unchanged)
    days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    day_of_year = zeros(size(month));
    
    for i = 1:length(month)
        if month(i) == 1
            day_of_year(i) = day(i);
        else
            day_of_year(i) = sum(days_in_month(1:month(i)-1)) + day(i);
        end
    end
end
function k_h_data_value = cal_k_h(z_height, surface_type)

    switch upper(surface_type)
    case 'A'
        kc = 1.174;a0 = 0.12;
    case 'B'
        kc = 1;a0 = 0.16;
    case 'C'
        kc = 0.785;a0 = 0.22;
    case 'D'
        kc = 0.564;a0 = 0.3;
    otherwise
        error('Surface type must be A, B, C or D');

    end
    k_h_data_value = kc*(z_height/10)^a0;
end

function gamma_t_value = interpolate_gamma_t(span_length, surface_type, filename)
% Interpolate to get gamma_t value
% Input parameters:
%   span_length: Bridge span (100~2000)
%   surface_type: Terrain roughness category ('A', 'B', 'C', 'D')
%   filename: Data file name (default: 'GAMMA_T.CSV')
% Output parameters:
%   gamma_t_value: Interpolated gamma_t value

    % Set default filename
    if nargin < 3
        filename = 'GAMMA_T.CSV';
    end
    
    % Read CSV file
    try
        data = readtable(filename);
    catch
        error('Unable to read file: %s. Please check if file exists.', filename);
    end
    
    % Check data validity
    if size(data, 2) < 5
        error('Insufficient columns in data file, at least 5 columns of data required');
    end
    
    % Extract span data
    spans = data{:, 1};
    
    % Check span range
    if span_length < min(spans) || span_length > max(spans)
        warning('Span %.1f exceeds data range [%.1f, %.1f], boundary values will be used', ...
                span_length, min(spans), max(spans));
    end
    
    % Select corresponding column based on terrain roughness category
    switch upper(surface_type)
        case 'A'
            gamma_t_data = data{:, 2};
        case 'B'
            gamma_t_data = data{:, 3};
        case 'C'
            gamma_t_data = data{:, 4};
        case 'D'
            gamma_t_data = data{:, 5};
        otherwise
            error('Surface type must be A, B, C or D');
    end
    
    % Perform interpolation
    gamma_t_value = interp1(spans, gamma_t_data, span_length, 'linear', 'extrap');
    
end

%% Drawing module
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