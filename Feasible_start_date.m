clear all; clc;
%% Input Parameters
% Erection duration (d)
d_tot = 100;   
%Bridge to be analysed
bridge_type = 1; % 1-the Severn Bridge; 2-the Xihoumen Bridge
% Cities to calculate
cal_series = [1, 2, 3, 4, 5, 6];% Input city numbers to calculate simultaneously
                            % 1:Dalian, 2:Qingdao, 3:Hangzhou, 4:Taipei, 5:Xiamen, 6:Hong Kong
%% Main function
switch bridge_type
    case 1
        span = 988; % main span of bridge (m)
        z = 40;     % reference height at bridge deck (m)
        global_params.U_f = [59.04,46.29,47.94,55.39,61.22,64.95,65.81,64.59,62.21,59.22]; % Flutter critical wind speeds for each stage
    case 2
        span = 1650;
        z = 50;
        global_params.U_f = [70.00,68.43,70.11,66.05,63.24,59.31,60.72,60.36,58.93,90.44]; % Flutter critical wind speeds for each stage
end

% Terrain roughness categories
surfaces = {'A', 'B', 'C', 'D'};
folder_path = 'WDSP_database';
city_names = {'Dalian', 'Qingdao', 'Hangzhou', 'Taipei', 'Xiamen', 'Hong Kong'};
file_names = {'wind_dalian_30.txt', 'wind_qingdao_30.txt', 'wind_hangzhou_30.txt', 'wind_taibei_30.txt', ...
    'wind_xiamen_30.txt', 'wind_hongkong_30.txt'};

% Parameter range
d_start_range = 1:1:365;   % Erection starting date (d)

% Initialize result matrices (24 rows × 365 columns)
feasibility_Results = zeros(24, length(d_start_range));
fail_count_Results = zeros(24, length(d_start_range));

% Progress display
total_iterations = length(cal_series) * length(surfaces) * length(d_start_range);
h = waitbar(0, 'Calculating...');
current_iteration = 0;

% Iterate through all cities and terrain roughness categories
for i_city = 1:length(cal_series)
    city_number = cal_series(i_city);
    
    if city_number >= 1 && city_number <= 6
        wind_file = file_names{city_number};
        city_name = city_names{city_number};
    else
        error('Invalid city number, please enter a number between 1-6');
    end
    
    % Read data
    data = read_wind_data(wind_file,folder_path);
    % Terrain correction
    U_corrected = terrain_correction(data, city_number);
    
    for i_surf = 1:length(surfaces)
        surface = surfaces{i_surf};
        
        % Calculate parameters
        k_f = 1.05;
        k_t = 1; % Terrain condition coefficient
        k_h = cal_k_h(z, surface);
        
        % Calculate Γ
        gamma_t = interpolate_gamma_t(span, surface, 'GAMMA_T.CSV');
        gamma_f = 1.15;  % Flutter stability partial factor
        gamma_a = 1.0;
        
        global_params.K = k_t * k_f * k_h;        % Coefficient K
        global_params.Gamma = gamma_t * gamma_f * gamma_a;    % Coefficient Gamma  
        global_params.N = length(global_params.U_f);          % Number of construction stages
        global_params.P_so = 0.9;
        global_params.P_se = global_params.P_so^(1/global_params.N);   
        global_params.city_name = city_name;
        
        % Initialize results for current city and terrain roughness category
        current_feasibility = zeros(1, length(d_start_range));
        current_fail_count = zeros(1, length(d_start_range));
        
        % Iterate through all starting dates
        for j = 1:length(d_start_range)
            d_start = d_start_range(j);
            
            % Check single scheme
            [result, fail_count] = check_single_scheme(data, U_corrected, d_tot, d_start, global_params);
            
            current_feasibility(j) = result;
            current_fail_count(j) = fail_count;
            
            % Update progress bar
            current_iteration = current_iteration + 1;
            if mod(current_iteration, 100) == 0
                progress = current_iteration / total_iterations;
                waitbar(progress, h, sprintf('Calculating... %.1f%%\nCity: %s, Surface: %s', progress*100, city_name, surface));
            end
        end
        
        % Calculate row position in result matrix (4 rows per city)
        result_row = (i_city - 1) * 4 + i_surf;
        feasibility_Results(result_row, :) = current_feasibility;
        fail_count_Results(result_row, :) = current_fail_count;
        fprintf('Completed: %s - Terrain roughness category %s\n', city_name, surface);
    end
end

close(h);
plot_all_feas_timeline(feasibility_Results, d_tot);

%% Sub functions
function [result, fail_count] = check_single_scheme(data,U_corrected,d_tot, d_start, global_params)
    % Function: Flutter safety check for a single scheme
    
    N = global_params.N;
    U_f = global_params.U_f;
    
    % Duration of each stage
    d = d_tot / N;
    
    fail_count = 0;
    
    % Check each stage
    for i = 1:N
        % Calculate stage start time
        d_s_i = d_start + (i-1) * d;
        % Normalize to 1-365
        d_s_i = mod(d_s_i - 1, 365) + 1;
        
        % Calculate flutter check wind speed for this stage
        U_c(i) = calculate_flutter_check_wind(data,U_corrected,d_s_i, d,global_params);
        
        % Compare with flutter critical wind speed
        if U_c(i) >= U_f(i)
            fail_count = fail_count + 1;
        end
    end
    
    % Determine result: pass only if all stages pass
    result = (fail_count == 0);
end

function U_c = calculate_flutter_check_wind(data,U_corrected,d_s, d,global_params)
    % Function 1: Calculate flutter check wind speed
    % Input: 
    %   d_s - stage start time (1~365)
    %   d - stage duration
    % Output: flutter check wind speed U_c
    years = unique(data.year);
    n_years = length(years); 
    
    % Get parameters for this stage
    [mu, alpha] = GM_params(data, U_corrected,d_s,d);
    
    % Extract global parameters
    K = global_params.K;
    Gamma = global_params.Gamma;
    P_se = global_params.P_se;
    
    % F(x) = exp(-exp(-(x-μ)/α))
    % F^(-1)(p) = μ - α * log(-log(p))
    y = - log(-log(P_se));
    K_c = sqrt(6)/pi*(y-0.5772);
    sigma = pi/sqrt(6)*alpha;
    F_inv = mu + alpha * y + sqrt(1+1.1*K_c^2+1.14*K_c)*sigma/sqrt(n_years);% Third term is finite sample correction
    
    % Calculate flutter check wind speed
    U_c = K * Gamma * F_inv;
end

function [mu, alpha] = GM_params(data, U_corrected, d_str, d_rp)
% Calculate Gumbel distribution parameters for the first stage
    EWS_history = group_first_period(data, U_corrected, d_str, d_rp);
    % Get all annual maximum wind speeds for the first stage (excluding 0 and NaN values)
    valid_data = EWS_history(EWS_history > 0 & ~isnan(EWS_history));
    
    if length(valid_data) < 2
        error('Insufficient valid data for the first stage, cannot calculate reliable parameters');
    end
    
    % Calculate mean and standard deviation
    mean_EWS = mean(valid_data);
    std_EWS = std(valid_data);
    
    % Calculate Gumbel parameters
    alpha = sqrt(6) / pi * std_EWS;
    mu = mean_EWS - 0.5772 * alpha;
end

function data = read_wind_data(filename,folder_path)
% Read wind speed data file
    full_filename = fullfile(folder_path, filename);
    
    fprintf('Reading data file: %s\n', full_filename);
    
    % Read text file
    raw_data = load(full_filename);
    
    % Extract each column data
    data.station = raw_data(:, 1);    % Station number
    data.year = raw_data(:, 2);       % Year
    data.month = raw_data(:, 3);      % Month
    data.day = raw_data(:, 4);        % Day
    data.V_600 = raw_data(:, 5)*0.5144;      % V_600 (10min daily maximum wind speed) original unit is knot
    data.V_3 = raw_data(:, 6)*0.5144;        % V_3 (3s daily maximum wind speed) original unit is knot
    
    fprintf('Data reading completed, total %d records\n', length(data.year));
end

% If beta needs to be calculated for each time point, it can be handled like this:
function U_corrected = terrain_correction(data,city_number)
% Terrain correction calculation

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
        
        % Filter valid data: V_3 and V_600 within reasonable ranges
        valid_idx = year_idx & (data.V_3 > 0) & (data.V_3 < 100) & ...
                   (data.V_600 > 5) & (data.V_600 < 50);
        
        % Calculate G values
        G_values = data.V_3(valid_idx) ./ data.V_600(valid_idx);
        
        % Check valid data days (excluding abnormal G values)
        valid_G = G_values(G_values > 0.1 & G_values < 2);
        
        if length(valid_G) >= 3  % At least 3 days of valid data per year
            yearly_G0(i) = mean(valid_G);
            valid_years(i) = true;
        else
            yearly_G0(i) = NaN;
            valid_years(i) = false;
        end
    end
    
    % Check if there are enough years for linear fitting
    valid_indices = find(valid_years);
    if length(valid_indices) < 3
        fprintf('Insufficient valid years (%d), no terrain correction applied\n', 3);
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
    

    % Apply terrain correction for each year
    for i = 1:n_years
        year_idx = (data.year == years(i));
        
        % Filter data that needs correction
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
    
    % Handle invalid data
    U_corrected(data.V_600 >= 50 | data.V_600 <= 0) = NaN;
    
    fprintf('Terrain correction completed, valid years %d/%d, fitting parameters: G(t)=%.4f*t+%.4f\n', ...
            length(valid_indices), n_years, p(1), p(2));
end

function EWS_history = group_first_period(data, U_corrected, d_str, d_rp)
% Calculate only the annual maximum wind speed for the first stage
    
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

        % Calculate start and end days of the first stage
        start_day = d_str;
        end_day = mod(d_str + d_rp - 1, days_in_year);
        
        % Determine the day range included in the stage
        if start_day <= end_day
            % No cross-year situation
            stage_idx = (year_day >= start_day) & (year_day <= end_day);
        else
            % Cross-year situation: from start_day to end of year, then from beginning of year to end_day
            stage_idx = (year_day >= start_day) | (year_day <= end_day);
        end
        
        % Extract wind speed data for this stage (excluding NaN values)
        stage_U = year_U(stage_idx);
        valid_stage_U = stage_U(~isnan(stage_U));
        
        if ~isempty(valid_stage_U)
            EWS_history(i) = max(valid_stage_U);
        else
            EWS_history(i) = NaN; % Set to NaN if no valid data
        end
    end
end

function day_of_year = calculate_day_of_year(month, day)
% Calculate the day of the year
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
    % Set default filename
    if nargin < 3
        filename = 'GAMMA_T.CSV';
    end
    
    % Read CSV file
    try
        data = readtable(filename);
    catch
        error('Cannot read file: %s. Please check if the file exists.', filename);
    end
    
    % Check data validity
    if size(data, 2) < 5
        error('Data file has insufficient columns, at least 5 columns of data are required');
    end
    
    % Extract span data
    spans = data{:, 1};
    
    % Check span range
    if span_length < min(spans) || span_length > max(spans)
        warning('Span %.1f is outside the data range [%.1f, %.1f], boundary values will be used', ...
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
function plot_all_feas_timeline(data,d_tot)
    % Correct data structure: 24 rows (6 cities × 4 categories) × 365 days
    cities = {'Dalian', 'Qingdao', 'Hangzhou', 'Taipei', 'Xiamen', 'Hong Kong'};
    categories = {'Category A', 'Category B', 'Category C', 'Category D'};
    
    % Create single plot
    figure('Position', [100, 100,950, 350]);
    
    % Set global font to Times
    set(gca, 'FontName', 'Times New Roman');
    
    % Define color scheme - use green shades for pass, red for fail
    pass_colors = [0.7, 1.0, 0.7;   % Light green - Category A
               0.6, 0.9, 0.6;   % Medium green - Category B
               0.3, 0.8, 0.3;   % Dark green - Category C
               0.1, 0.6, 0.1];  % Very dark green - Category D
    
    fail_color = [0.9, 0.2, 0.2];   % Red - Fail
    
    hold on;
    
    % Plot data for each city and category
    for city = 1:6
    for category = 1:4
        % Calculate data row index: each city occupies 4 rows, arranged by category order
        row_index = (city-1)*4 + category;
        city_category_data = data(row_index, :);
        
        % Calculate y coordinate: each city occupies a height range, Category A~D from top to bottom
        y_base = (6 - city) * 5;  % Dalian at top (y=25), Hong Kong at bottom (y=5)
        y_pos = y_base - category + 1;  % Category A at top, D at bottom
        
        % Find pass and fail days
        pass_days = find(city_category_data == 1);
        fail_days = find(city_category_data == 0);
        
        % Plot pass points (using different shades of green)
        if ~isempty(pass_days)
            h_pass(category) = scatter(pass_days, y_pos * ones(size(pass_days)), ...
                    20, pass_colors(category, :), 'filled', ...
                    'Marker', 'o', 'MarkerEdgeColor', 'none', ...
                    'DisplayName', categories{category});
        end
        
        % Plot fail points (red solid circles) - only add to legend when first plotted
        if ~isempty(fail_days)
            if city == 1 && category == 1
                h_fail = scatter(fail_days, y_pos * ones(size(fail_days)), ...
                        20, fail_color, 'filled', ...
                        'Marker', 'o', 'MarkerEdgeColor', 'none', ...
                        'DisplayName', 'Fail');
            else
                scatter(fail_days, y_pos * ones(size(fail_days)), ...
                        20, fail_color, 'filled', ...
                        'Marker', 'o', 'MarkerEdgeColor', 'none', ...
                        'HandleVisibility', 'off');
            end
        end
    end
    end
    
    % Set axes and labels
    ylim([-4.5, 26.5]);
    xlim([0, 368]);
    xlabel('$d_{str} (d)$', 'FontSize', 15, 'FontWeight', 'normal', 'FontName', 'Times New Roman','Interpreter','latex');
    box on
    % Add city separation lines and labels
    for city = 1:6
    y_city = (6 - city) * 5 + 0.5;
    % City separation line
    plot([1, 365], [y_city, y_city], 'k-', 'LineWidth', 0.5, 'Color', [0.5, 0.5, 0.5]);
    
    % City name labels
    text(-45, y_city - 2, cities{city}, 'FontSize', 12, 'FontWeight', 'normal', ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
         'FontName', 'Times New Roman');
    end
    % Add d_tot annotation at top right
text(-47, 27.5, sprintf('($d_{tot} $= %dd)', d_tot), 'FontSize', 14, 'FontWeight', 'normal', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
     'FontName', 'Times New Roman','Interpreter','latex');
    
    % Add grid
    grid on;
    set(gca, 'GridAlpha', 0.2, 'XGrid', 'on', 'YGrid', 'off');
    set(gcf, 'Color', 'white');
    set(gca, 'FontSize', 12, 'FontName', 'Times');
    % Set coordinate axis line width to 1.5 points
    set(gca, 'LineWidth', 1.5);
    % Set y-axis ticks
    yticks([]);
    xticks(1:60:365);  % One tick every 60 days
    xticklabels({'1', '60', '120', '180', '240', '300', '360'});
    
    % Create legend (only 5 types: 4 Category Pass + 1 Fail)
    legend([h_pass(1), h_pass(2), h_pass(3), h_pass(4), h_fail], ...
       {'Category A Safe', 'Category B Safe', 'Category C Safe', 'Category D Safe', 'Unsafe'}, ...
       'Location', 'northeast', 'FontSize', 12, 'FontName', 'Times New Roman',NumColumns=5,Box='off');
    
    hold off;
end