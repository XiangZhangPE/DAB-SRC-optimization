% clc, close all;
%% Define constants
fs = 250e3; % switching frequency, user defined
Lr = 12e-6; % resonant inductor, user defined
Cr = 48e-9; % resonant capacitor, user defined
Z0 = sqrt(Lr/Cr);
fr = 1/(2*pi*sqrt(Lr*Cr));
r = fr/fs;

Vin = 100; % input voltage, user defined
Vo = 70;  % output voltage, user defined

% unused variables
omega_r = 2*pi*fr;
td = 50e-9;
Dd = td*fs;
Ln = 100;
ILm0 = pi*r*Vo/(2 * Ln * Z0);
Ioss = 0.06;

% % plot background ZVS figure
% handle = figure(3);
% plotoptions = 2;
% ZVSPlot_DABsrc(r,Vo,Vin,Z0, plotoptions); hold on;

%% Optimization problem
DoF = 4;

% Define bounds for optimization variables based on DoF
if DoF == 1
    bounds = [0; 0.5];  % Lower and upper bounds for Dp
elseif DoF == 2
    bounds = [0 0; 0.5 1];  % Bounds for Dp, Dy1
elseif DoF == 3
    bounds = [0 0 0; 0.5 1 1];  % Bounds for Dp, Dy1, Dy2
elseif DoF == 4
    bounds = [0 0 0 1/1.3; 0.5 1 1 1/1.05];  % Bounds for Dp, Dy1, Dy2, r
end

% Update optimization options to increase iterations
options = optimoptions('fmincon', ...
    'Display', 'off', ...
    'MaxIterations', 2000, ...  % Increase maximum iterations
    'MaxFunctionEvaluations', 10000, ...  % Increase function evaluations
    'OptimalityTolerance', 2e-9, ...
    'StepTolerance', 1e-10);

% Define the power ranges and output votage ranges
low_power_range = linspace(0.025, 0.45, 25);
high_power_range = linspace(0.45, 1.2, 25);
PoN_values = [low_power_range, high_power_range]; % Include both power ranges
Vo_values = 60:2.5:135; % Adjust the range as needed

% Define different initial conditions for each region
% low_power_buck_initial = [0.012, 0.29, 0.38, 0.77];
% high_power_buck_initial = [0.19, 0.96, 1, 0.858];
% low_power_boost_initial = [0.012, 0.38, 0.29, 0.77];
% high_power_boost_initial = [0.175, 1, 1, 0.868];
% low_power_unity_initial = [0.003, 1, 1, 0.77];
% high_power_unity_initial = [0.07, 1, 1, 0.95];

low_power_buck_initial = [0.0259, 0.2315, 0.4178, 0.7692];
high_power_buck_initial =[0.2665, 0.7360, 1, 0.8349];
low_power_boost_initial =[0.0187, 0.4456, 0.3110, 0.7692];
high_power_boost_initial=[0.2166, 1, 0.8764, 0.8477];
low_power_unity_initial =[0.003, 1, 1, 0.77];
high_power_unity_initial=[0.07, 1, 1, 0.95];

% 定义欧几里得距离的阈值
distance_threshold = 0.02; % 两个点之间的最大允许距离

% Preallocate 3D arrays for storing results
optimal_points = zeros(length(PoN_values), DoF, length(Vo_values));
optimal_power = zeros(length(PoN_values), length(Vo_values));
optimal_Irms = zeros(length(PoN_values), length(Vo_values));
Ioffp12 = zeros(length(PoN_values), length(Vo_values));
Ioffp34 = zeros(length(PoN_values), length(Vo_values));
Ioffs12 = zeros(length(PoN_values), length(Vo_values));
Ioffs34 = zeros(length(PoN_values), length(Vo_values));
ILmax = zeros(length(PoN_values), length(Vo_values));
actual_PoN_values = zeros(length(PoN_values), length(Vo_values));

% Use parfor to parallelize the Vo loop
parfor v_idx = 1:length(Vo_values)
    Vo = Vo_values(v_idx); % Update Vo for this iteration

    % Local variables to store results for this Vo
    local_optimal_points = zeros(length(PoN_values), DoF);
    local_optimal_power = zeros(length(PoN_values), 1);
    local_optimal_Irms = zeros(length(PoN_values), 1);
    local_Ioffp12 = zeros(length(PoN_values), 1);
    local_Ioffp34 = zeros(length(PoN_values), 1);
    local_Ioffs12 = zeros(length(PoN_values), 1);
    local_Ioffs34 = zeros(length(PoN_values), 1);
    local_ILmax = zeros(length(PoN_values), 1);
    local_actual_PoN_values = zeros(length(PoN_values), 1);

    % Determine initial conditions based on Vo
    if (Vo > 1.05 * Vin)
        high_power_initial = high_power_boost_initial;
        low_power_initial = low_power_boost_initial;
    elseif (Vo < 0.95 * Vin)
        high_power_initial = high_power_buck_initial;
        low_power_initial = low_power_buck_initial;
    else
        high_power_initial = high_power_unity_initial;
        low_power_initial = low_power_unity_initial;
    end

    % For each power target
    for i = 1:length(PoN_values)
        PoN_target = PoN_values(i);

        % Define objective function with penalty
        obj_fun = @(vars) J_function(DoF, vars, Vin, Vo, PoN_target);

        % Define nonlinear constraints
        nonlcon = @(vars) zvs_function(DoF, vars, Vin, Vo, Z0);

        % Perform separate optimizations for low and high power initial guesses
        [optimal_vars_low, fval_low, exitflag_low, ~] = fmincon(obj_fun, low_power_initial, [], [], [], [], bounds(1,:), bounds(2,:), nonlcon, options);
        [optimal_vars_high, fval_high, exitflag_high, ~] = fmincon(obj_fun, high_power_initial, [], [], [], [], bounds(1,:), bounds(2,:), nonlcon, options);

        % Select the better result
        if exitflag_low > 0 && exitflag_high > 0
            if fval_low <= fval_high
                optimal_vars = optimal_vars_low;
            else
                optimal_vars = optimal_vars_high;
            end
        elseif exitflag_low > 0
            optimal_vars = optimal_vars_low;
        elseif exitflag_high > 0
            optimal_vars = optimal_vars_high;
        else
            warning('⚠️ Both optimizations did NOT converge for Vo = %.2f, PoN_target = %.4f \n', Vo, PoN_target);
            continue; % Skip to next iteration
        end

        % Store results in local variables
        local_optimal_points(i, :) = optimal_vars;
        local_optimal_power(i) = f_function(DoF, optimal_vars, Vin, Vo);
        local_optimal_Irms(i) = g_function(DoF, optimal_vars, Vin, Vo);
        [Ip12_val, Ip34_val, Is12_val, Is34_val, ILmax_val] = h_function(DoF, optimal_vars, Vin, Vo);
        local_Ioffp12(i) = Ip34_val;
        local_Ioffp34(i) = Ip12_val;
        local_Ioffs12(i) = Is34_val;
        local_Ioffs34(i) = Is12_val;
        local_ILmax(i) = ILmax_val;
        local_actual_PoN_values(i) = f_function(DoF, optimal_vars, Vin, Vo);

        % Update initial guesses for the next iteration
        low_power_initial = optimal_vars;
        high_power_initial = optimal_vars;
    end

    % Assign local results to global arrays
    optimal_points(:, :, v_idx) = local_optimal_points;
    optimal_power(:, v_idx) = local_optimal_power;
    optimal_Irms(:, v_idx) = local_optimal_Irms;
    Ioffp12(:, v_idx) = local_Ioffp12;
    Ioffp34(:, v_idx) = local_Ioffp34;
    Ioffs12(:, v_idx) = local_Ioffs12;
    Ioffs34(:, v_idx) = local_Ioffs34;
    ILmax(:, v_idx) = local_ILmax;
    actual_PoN_values(:, v_idx) = local_actual_PoN_values;
end

% parfor循环完成后调用
[optimal_points, optimal_power, optimal_Irms, Ioffp12, Ioffp34, Ioffs12, Ioffs34, ILmax, actual_PoN_values] = ...
    remove_outlier_points(optimal_points, optimal_power, optimal_Irms, Ioffp12, Ioffp34, Ioffs12, Ioffs34, ILmax, actual_PoN_values, Vin, Vo_values, DoF);

% 然后生成网格和绘图
[PoN_grid, Vo_grid] = meshgrid(PoN_values, Vo_values / Vin);

% Figure 1: Normalized turn-off current
figure(1);
hold on;
mesh(PoN_grid, Vo_grid, Ioffp12', 'FaceColor', 'r', 'EdgeColor', [0.3 0.3 0.3], 'EdgeAlpha', 0.3,  'FaceAlpha', 0.4, 'DisplayName', 'I_{off,P12}');
mesh(PoN_grid, Vo_grid, Ioffp34', 'FaceColor', 'b', 'EdgeColor', [0.3 0.3 0.3], 'EdgeAlpha', 0.3,  'FaceAlpha', 0.4, 'DisplayName', 'I_{off,P34}');
mesh(PoN_grid, Vo_grid, Ioffs12', 'FaceColor', 'm', 'EdgeColor', [0.3 0.3 0.3], 'EdgeAlpha', 0.3,  'FaceAlpha', 0.4, 'DisplayName', 'I_{off,S12}');
mesh(PoN_grid, Vo_grid, Ioffs34', 'FaceColor', 'c', 'EdgeColor', [0.3 0.3 0.3], 'EdgeAlpha', 0.3,  'FaceAlpha', 0.2, 'DisplayName', 'I_{off,S34}');
mesh(PoN_grid, Vo_grid, ILmax', 'FaceColor', 'y', 'EdgeColor', [0.3 0.3 0.3], 'EdgeAlpha', 0.1,  'DisplayName', 'I_{L,max}');
xlabel('Normalized power PoN = Pout/PN', 'Interpreter', 'latex');
ylabel('Voltage ratio M = Vo/Vin', 'Interpreter', 'latex');
zlabel('Normalized turn-off current', 'Interpreter', 'latex');
title('Ioff under different Power and Vo', 'Interpreter', 'latex');
legend('show');
grid on;
view(3); % 3D view

% Figure 2: Phase-shift ratios
figure(2);
hold on;

% Subplot 1: Dp
subplot(2, 2, 1);
mesh(PoN_grid, Vo_grid, squeeze(optimal_points(:, 1, :))', 'FaceColor', 'r', 'EdgeColor', [0.3 0.3 0.3], 'EdgeAlpha', 0.3,  'FaceAlpha', 0.5);
xlabel('Normalized power PoN = Pout/PN', 'Interpreter', 'latex');
ylabel('Voltage ratio M = Vo/Vin', 'Interpreter', 'latex');
zlabel('Dp', 'Interpreter', 'latex');
title('Outer-phase-shift ratio Dp', 'Interpreter', 'latex');
grid on;
view(3);

% Subplot 2: Dy1
subplot(2, 2, 2);
mesh(PoN_grid, Vo_grid, squeeze(optimal_points(:, 2, :))', 'FaceColor', 'g', 'EdgeColor', [0.3 0.3 0.3], 'EdgeAlpha', 0.3,  'FaceAlpha', 0.5);
xlabel('Normalized power PoN = Pout/PN', 'Interpreter', 'latex');
ylabel('Voltage ratio M = Vo/Vin', 'Interpreter', 'latex');
zlabel('Dy1', 'Interpreter', 'latex');
title('Primary Inner-phase-shift duty Dy1', 'Interpreter', 'latex');
grid on;
view(3);

% Subplot 3: Dy2
subplot(2, 2, 3);
mesh(PoN_grid, Vo_grid, squeeze(optimal_points(:, 3, :))', 'FaceColor', 'b', 'EdgeColor', [0.3 0.3 0.3], 'EdgeAlpha', 0.3,  'FaceAlpha', 0.5);
xlabel('Normalized power PoN = Pout/PN', 'Interpreter', 'latex');
ylabel('Voltage ratio M = Vo/Vin', 'Interpreter', 'latex');
zlabel('Dy2', 'Interpreter', 'latex');
title('Secondary Inner-phase-shift duty Dy2', 'Interpreter', 'latex');
grid on;
view(3);

% Subplot 4: r (only if DoF >= 4)
if DoF >= 4
    subplot(2, 2, 4);
    mesh(PoN_grid, Vo_grid, squeeze(optimal_points(:, 4, :))', 'FaceColor', 'y', 'EdgeColor', [0.3 0.3 0.3], 'EdgeAlpha', 0.3,  'FaceAlpha', 0.5);
    xlabel('Normalized power PoN = Pout/PN', 'Interpreter', 'latex');
    ylabel('Voltage ratio M = Vo/Vin', 'Interpreter', 'latex');
    zlabel('r', 'Interpreter', 'latex');
    title('Frequency ratio r = fr/fs', 'Interpreter', 'latex');
    grid on;
    view(3);
end

sgtitle('Optimized Phase-Shift Ratios', 'Interpreter', 'latex');






% J_function
function J = J_function(DoF, vars, Vin, Vo, PoN_target)
    % calculate PoN
    PoN = f_function(DoF, vars, Vin, Vo);
    % calculate Irms
    Irms = g_function(DoF, vars, Vin, Vo);
    % calculate Ioff (you can ignore the return value if not needed)
    h_function(DoF, vars, Vin, Vo);
    
    % construct cost function
    lambda = 800;  % add penalty coefficient (weight λ adjustable)
    penalty = lambda * (PoN - PoN_target)^2;
    
    % final cost function
    J = Irms + penalty;
end


function [Dp, Dy1, Dy2, r] = assign_variables(DoF, vars, Vin, Vo)
    % Default values
    r = 1/1.2;
    
    % Assign based on DoF
    if DoF >= 1
        Dp = vars(1); Dy1 = 1; Dy2 = 1;
    end
    
    if DoF >= 2
        if (Vin > Vo)
            Dy1 = vars(2);
            Dy2 = 1;
        else
            Dy1 = 1;
            Dy2 = vars(2);
        end
    end
    
    if DoF >= 3
        Dy1 = vars(2); Dy2 = vars(3);
    end
    
    if DoF >= 4
        r = vars(4);
    end
end


% Define power function
function f = f_function(DoF, vars, Vin, Vo)

    % Get variables based on DoF
    [Dp, Dy1, Dy2, r] = assign_variables(DoF, vars, Vin, Vo);
    
    region_flag((Dy1>Dy2) & (Dp<(Dy1-Dy2)/2)) = 1;  %  A
    region_flag((Dy2>=Dy1) & (Dp<(Dy2-Dy1)/2)) = 2;  %  B
    region_flag((Dp>=abs(Dy1-Dy2)/2) & (Dp<(Dy2+Dy1)/2) & (Dp<(1-(Dy2+Dy1)/2))) = 3;  %  C
    region_flag((((Dy2+Dy1)/2)<= Dp) & (Dp<(1-(Dy2+Dy1)/2))) = 4;  %  D
    region_flag(((1-(Dy2+Dy1)/2)<= Dp) & (((Dy2+Dy1)/2)>= Dp)) = 5;  %  E
    % region_flag(((1-(Dy2-Dy1)/2)<=Dp) & (Dy2>Dy1)) = 6;  %  F
    % region_flag(((1-(Dy1-Dy2)/2)<=Dp) & (Dy2<=Dy1)) = 7;  %  G
    % region_flag((Dp>=(1-(Dy1+Dy2)/2)) & (Dp>=(Dy1+Dy2)/2) & (Dp<(1-abs(Dy1-Dy2)/2))) = 8;  %  H

    switch region_flag
        case 1
            % Region A
            f =  1/2.*((sin(pi*r.*Dy2/2).*sin(pi*r.*Dp).*cos(pi*r.*(1-Dy1)/2))./(r.*cos(pi*r/2)));
        case 2
            % Region B
            f = 1/2.*((sin(pi*r.*Dy1/2).*sin(pi*r.*Dp).*cos(pi*r.*(1-Dy2)/2))./(r.*cos(pi*r/2)));
        case 3
            % Region C
            f = 1/(4).*(((cos(pi.*r.*(Dy2-Dy1)/2).*cos(pi.*r/2.*(2.*Dp-1)) + sin(pi.*r.*(Dy1+Dy2-1)/2).*sin(pi.*r.*(Dp)))./(r.*cos(pi.*r./2)) - 1./r));
        case 4
            % Region D
            f = 1/2.*(((sin(pi*r.*Dy1/2).*cos(pi*r.*(2*Dp-1)/2).*sin(pi*r.*Dy2/2))./(r.*cos(pi*r/2))));
        case 5
            % Region E
            f = 1/2.*((cos(pi.*r.*(Dp-0.5)).*cos(pi.*r.*(1-Dy1)./2).*cos(pi.*r.*(1-Dy2)./2))./(r.*cos(pi.*r./2))-1./(r));
        % case 6
        %     % Region F
        %     f = 1/2.*((sin(pi*r.*Dy1/2).*sin(pi*r.*(1-Dp)).*cos(pi*r.*(1-Dy2)/2))./(r.*cos(pi*r/2)));
        % case 7
        %     % Region G
        %     f = 1/2.*((sin(pi*r.*Dy2/2).*sin(pi*r.*(1-Dp)).*cos(pi*r.*(1-Dy1)/2))./(r.*cos(pi*r/2)));
        % case 8
        %     % Region H
        %     f = 1/(4).*(((cos(pi.*r.*(Dy2-Dy1)/2).*cos(pi.*r/2.*(2.*Dp-1)) + sin(pi.*r.*(Dy1+Dy2-1)/2).*sin(pi.*r.*(1-Dp)))./(r.*cos(pi.*r./2)) - 1./r));
        otherwise
            error('Unknown region');
    end

end



% Define Irms function
function g = g_function(DoF, vars, Vin, Vo)
    % Get variables based on DoF
    [Dp, Dy1, Dy2, r] = assign_variables(DoF, vars, Vin, Vo);

    region_flag((Dy1>Dy2) & (Dp<(Dy1-Dy2)/2)) = 1;  %  A
    region_flag((Dy2>=Dy1) & (Dp<(Dy2-Dy1)/2)) = 2;  %  B
    region_flag((Dp>=abs(Dy1-Dy2)/2) & (Dp<(Dy2+Dy1)/2) & (Dp<(1-(Dy2+Dy1)/2))) = 3;  %  C
    region_flag((((Dy2+Dy1)/2)<= Dp) & (Dp<(1-(Dy2+Dy1)/2))) = 4;  %  D
    region_flag(((1-(Dy2+Dy1)/2)<= Dp) & (((Dy2+Dy1)/2)>= Dp)) = 5;  %  E
    % region_flag(((1-(Dy2-Dy1)/2)<=Dp) & (Dy2>Dy1)) = 6;  %  F
    % region_flag(((1-(Dy1-Dy2)/2)<=Dp) & (Dy2<=Dy1)) = 7;  %  G
    % region_flag((Dp>=(1-(Dy1+Dy2)/2)) & (Dp>=(Dy1+Dy2)/2) & (Dp<(1-abs(Dy1-Dy2)/2))) = 8;  %  H

    % calulate objective function
    switch region_flag
        case 1
            % Region A
            g = 1./Vin.*sqrt((1/2/pi./r).*((2.*Vin.^2.*sin(pi.*r).*cos((Dy1.*pi.*r)./2).^2-Vo.^2.*sin(Dy2.*pi.*r)-2.*Vin.^2.*sin(pi.*r)-2.*Vo.^2.*sin(pi.*r)-Vin.^2.*sin(Dy1.*pi.*r)+2.*Vo.^2.*sin(pi.*r).*cos((Dy2.*pi.*r)./2).^2+2.*Vin.^2.*pi.*r+2.*Vo.^2.*pi.*r-Dy1.*Vin.^2.*pi.*r-Dy2.*Vo.^2.*pi.*r-2.*Vin.^2.*pi.*r.*cos((Dy1.*pi.*r)./2).^2-2.*Vo.^2.*pi.*r.*cos((Dy2.*pi.*r)./2).^2-2.*Vin.^2.*cos(pi.*r).*cos((Dy1.*pi.*r)./2).*sin((Dy1.*pi.*r)./2)-2.*Vo.^2.*cos(pi.*r).*cos((Dy2.*pi.*r)./2).*sin((Dy2.*pi.*r)./2)+2.*Dy1.*Vin.^2.*pi.*r.*cos((Dy1.*pi.*r)./2).^2+2.*Dy2.*Vo.^2.*pi.*r.*cos((Dy2.*pi.*r)./2).^2+4.*Vin.*Vo.*cos(Dp.*pi.*r).*cos((Dy1.*pi.*r)./2).*sin((Dy2.*pi.*r)./2)-Dy1.*Vin.^2.*pi.*r.*cos(pi.*r)-Dy2.*Vo.^2.*pi.*r.*cos(pi.*r)+2.*Dy1.*Vin.^2.*pi.*r.*cos(pi.*r).*cos((Dy1.*pi.*r)./2).^2+2.*Dy2.*Vo.^2.*pi.*r.*cos(pi.*r).*cos((Dy2.*pi.*r)./2).^2+4.*Vin.*Vo.*cos(pi.*r).*cos(Dp.*pi.*r).*cos((Dy1.*pi.*r)./2).*sin((Dy2.*pi.*r)./2)+4.*Vin.*Vo.*sin(pi.*r).*cos(Dp.*pi.*r).*sin((Dy1.*pi.*r)./2).*sin((Dy2.*pi.*r)./2)-4.*Vin.*Vo.*pi.*r.*cos(Dp.*pi.*r).*sin((Dy1.*pi.*r)./2).*sin((Dy2.*pi.*r)./2)+2.*Dy1.*Vin.^2.*pi.*r.*sin(pi.*r).*cos((Dy1.*pi.*r)./2).*sin((Dy1.*pi.*r)./2)+2.*Dy2.*Vo.^2.*pi.*r.*sin(pi.*r).*cos((Dy2.*pi.*r)./2).*sin((Dy2.*pi.*r)./2)+4.*Dp.*Vin.*Vo.*pi.*r.*cos((Dy1.*pi.*r)./2).*sin(Dp.*pi.*r).*sin((Dy2.*pi.*r)./2)+2.*Dy1.*Vin.*Vo.*pi.*r.*cos(Dp.*pi.*r).*sin((Dy1.*pi.*r)./2).*sin((Dy2.*pi.*r)./2)-2.*Dy2.*Vin.*Vo.*pi.*r.*cos(Dp.*pi.*r).*cos((Dy1.*pi.*r)./2).*cos((Dy2.*pi.*r)./2)-2.*Dy2.*Vin.*Vo.*pi.*r.*cos(pi.*r).*cos(Dp.*pi.*r).*cos((Dy1.*pi.*r)./2).*cos((Dy2.*pi.*r)./2)+4.*Dp.*Vin.*Vo.*pi.*r.*cos(pi.*r).*cos((Dy1.*pi.*r)./2).*sin(Dp.*pi.*r).*sin((Dy2.*pi.*r)./2)+2.*Dy1.*Vin.*Vo.*pi.*r.*cos(pi.*r).*cos(Dp.*pi.*r).*sin((Dy1.*pi.*r)./2).*sin((Dy2.*pi.*r)./2)-2.*Dy1.*Vin.*Vo.*pi.*r.*sin(pi.*r).*cos(Dp.*pi.*r).*cos((Dy1.*pi.*r)./2).*sin((Dy2.*pi.*r)./2)-2.*Dy2.*Vin.*Vo.*pi.*r.*sin(pi.*r).*cos(Dp.*pi.*r).*cos((Dy2.*pi.*r)./2).*sin((Dy1.*pi.*r)./2)+4.*Dp.*Vin.*Vo.*pi.*r.*sin(pi.*r).*sin(Dp.*pi.*r).*sin((Dy1.*pi.*r)./2).*sin((Dy2.*pi.*r)./2))./(2.*1.^2.*(cos(pi.*r)+1))));
        case 2
            % Region B
            g = 1./Vin.*sqrt((1/2/pi./r).*((2.*Vin.^2.*sin(pi.*r).*cos((Dy1.*pi.*r)./2).^2 - Vo.^2.*sin(Dy2.*pi.*r) - 2.*Vin.^2.*sin(pi.*r) - 2.*Vo.^2.*sin(pi.*r) - Vin.^2.*sin(Dy1.*pi.*r) + 2.*Vo.^2.*sin(pi.*r).*cos((Dy2.*pi.*r)./2).^2 + 2.*Vin.^2.*pi.*r + 2.*Vo.^2.*pi.*r - Dy1.*Vin.^2.*pi.*r - Dy2.*Vo.^2.*pi.*r - 2.*Vin.^2.*pi.*r.*cos((Dy1.*pi.*r)./2).^2 - 2.*Vo.^2.*pi.*r.*cos((Dy2.*pi.*r)./2).^2 - 2.*Vin.^2.*cos(pi.*r).*cos((Dy1.*pi.*r)./2).*sin((Dy1.*pi.*r)./2) - 2.*Vo.^2.*cos(pi.*r).*cos((Dy2.*pi.*r)./2).*sin((Dy2.*pi.*r)./2) + 2.*Dy1.*Vin.^2.*pi.*r.*cos((Dy1.*pi.*r)./2).^2 + 2.*Dy2.*Vo.^2.*pi.*r.*cos((Dy2.*pi.*r)./2).^2 + 4.*Vin.*Vo.*cos(Dp.*pi.*r).*cos((Dy2.*pi.*r)./2).*sin((Dy1.*pi.*r)./2) - Dy1.*Vin.^2.*pi.*r.*cos(pi.*r) - Dy2.*Vo.^2.*pi.*r.*cos(pi.*r) + 2.*Dy1.*Vin.^2.*pi.*r.*cos(pi.*r).*cos((Dy1.*pi.*r)./2).^2 + 2.*Dy2.*Vo.^2.*pi.*r.*cos(pi.*r).*cos((Dy2.*pi.*r)./2).^2 + 4.*Vin.*Vo.*cos(pi.*r).*cos(Dp.*pi.*r).*cos((Dy2.*pi.*r)./2).*sin((Dy1.*pi.*r)./2) + 4.*Vin.*Vo.*sin(pi.*r).*cos(Dp.*pi.*r).*sin((Dy1.*pi.*r)./2).*sin((Dy2.*pi.*r)./2) - 4.*Vin.*Vo.*pi.*r.*cos(Dp.*pi.*r).*sin((Dy1.*pi.*r)./2).*sin((Dy2.*pi.*r)./2) + 2.*Dy1.*Vin.^2.*pi.*r.*sin(pi.*r).*cos((Dy1.*pi.*r)./2).*sin((Dy1.*pi.*r)./2) + 2.*Dy2.*Vo.^2.*pi.*r.*sin(pi.*r).*cos((Dy2.*pi.*r)./2).*sin((Dy2.*pi.*r)./2) + 4.*Dp.*Vin.*Vo.*pi.*r.*cos((Dy2.*pi.*r)./2).*sin(Dp.*pi.*r).*sin((Dy1.*pi.*r)./2) + 2.*Dy2.*Vin.*Vo.*pi.*r.*cos(Dp.*pi.*r).*sin((Dy1.*pi.*r)./2).*sin((Dy2.*pi.*r)./2) - 2.*Dy1.*Vin.*Vo.*pi.*r.*cos(Dp.*pi.*r).*cos((Dy1.*pi.*r)./2).*cos((Dy2.*pi.*r)./2) - 2.*Dy1.*Vin.*Vo.*pi.*r.*cos(pi.*r).*cos(Dp.*pi.*r).*cos((Dy1.*pi.*r)./2).*cos((Dy2.*pi.*r)./2) + 4.*Dp.*Vin.*Vo.*pi.*r.*cos(pi.*r).*cos((Dy2.*pi.*r)./2).*sin(Dp.*pi.*r).*sin((Dy1.*pi.*r)./2) - 2.*Dy1.*Vin.*Vo.*pi.*r.*sin(pi.*r).*cos(Dp.*pi.*r).*cos((Dy1.*pi.*r)./2).*sin((Dy2.*pi.*r)./2) + 2.*Dy2.*Vin.*Vo.*pi.*r.*cos(pi.*r).*cos(Dp.*pi.*r).*sin((Dy1.*pi.*r)./2).*sin((Dy2.*pi.*r)./2) - 2.*Dy2.*Vin.*Vo.*pi.*r.*sin(pi.*r).*cos(Dp.*pi.*r).*cos((Dy2.*pi.*r)./2).*sin((Dy1.*pi.*r)./2) + 4.*Dp.*Vin.*Vo.*pi.*r.*sin(pi.*r).*sin(Dp.*pi.*r).*sin((Dy1.*pi.*r)./2).*sin((Dy2.*pi.*r)./2))./(2.*1.^2.*(cos(pi.*r) + 1))));

        case 3
            % Region C
            g = 0;
        case 4
            % Region D
            g = 0;
        case 5
            % Region E
            g = 1./Vin.*sqrt((1/2/pi./r).*(-(2.*Vin.^2.*sin(Dy1.*pi.*r) + 2.*Vo.^2.*sin(Dy2.*pi.*r) + 2.*Vin.^2.*sin(pi.*r.*(Dy1 - 1)) + 2.*Vo.^2.*sin(pi.*r.*(Dy2 - 1)) + 2.*Vin.^2.*sin(pi.*r) + 2.*Vo.^2.*sin(pi.*r) + 2.*Vin.*Vo.*sin((pi.*r.*(2.*Dp - Dy1 - Dy2 + 2))./2) - 2.*Vin.*Vo.*sin((pi.*r.*(Dy1 - 2.*Dp + Dy2))./2) + 2.*Vin.*Vo.*sin((pi.*r.*(2.*Dp + Dy1 + Dy2 - 2))./2) + 2.*Vin.*Vo.*sin((pi.*r.*(2.*Dp + Dy1 + Dy2 - 4))./2) + 2.*Vin.*Vo.*sin((pi.*r.*(2.*Dp + Dy1 - Dy2))./2) + 2.*Vin.*Vo.*sin((pi.*r.*(2.*Dp - Dy1 + Dy2))./2) + 2.*Vin.*Vo.*sin((pi.*r.*(2.*Dp + Dy1 - Dy2 - 2))./2) + 2.*Vin.*Vo.*sin((pi.*r.*(2.*Dp - Dy1 + Dy2 - 2))./2) - 2.*Vin.^2.*pi.*r - 2.*Vo.^2.*pi.*r + 2.*Vin.^2.*pi.*r.*cos(Dy1.*pi.*r) + 2.*Vo.^2.*pi.*r.*cos(Dy2.*pi.*r) - 2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(Dy1 - 2.*Dp + Dy2))./2) + 4.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp + Dy1 + Dy2 - 2))./2) + 2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp + Dy1 + Dy2 - 4))./2) - 2.*Dy1.*Vin.^2.*pi.*r.*cos(Dy1.*pi.*r) - 2.*Dy2.*Vo.^2.*pi.*r.*cos(Dy2.*pi.*r) + 2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp + Dy1 - Dy2))./2) + 2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp - Dy1 + Dy2))./2) - 2.*Dy1.*Vin.^2.*pi.*r.*cos(pi.*r.*(Dy1 - 1)) - 2.*Dy2.*Vo.^2.*pi.*r.*cos(pi.*r.*(Dy2 - 1)) - 2.*Dp.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp - Dy1 - Dy2 + 2))./2) + Dy1.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp - Dy1 - Dy2 + 2))./2) + Dy2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp - Dy1 - Dy2 + 2))./2) - 2.*Dp.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(Dy1 - 2.*Dp + Dy2))./2) - 2.*Dp.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp + Dy1 + Dy2 - 2))./2) - 2.*Dp.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp + Dy1 + Dy2 - 4))./2) + Dy1.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(Dy1 - 2.*Dp + Dy2))./2) + Dy2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(Dy1 - 2.*Dp + Dy2))./2) - Dy1.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp + Dy1 + Dy2 - 2))./2) - Dy2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp + Dy1 + Dy2 - 2))./2) - Dy1.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp + Dy1 + Dy2 - 4))./2) - Dy2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp + Dy1 + Dy2 - 4))./2) - 2.*Dp.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp + Dy1 - Dy2))./2) - 2.*Dp.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp - Dy1 + Dy2))./2) - 2.*Dp.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp + Dy1 - Dy2 - 2))./2) - 2.*Dp.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp - Dy1 + Dy2 - 2))./2) - Dy1.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp + Dy1 - Dy2))./2) + Dy1.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp - Dy1 + Dy2))./2) + Dy2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp + Dy1 - Dy2))./2) - Dy2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp - Dy1 + Dy2))./2) - Dy1.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp + Dy1 - Dy2 - 2))./2) + Dy1.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp - Dy1 + Dy2 - 2))./2) + Dy2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp + Dy1 - Dy2 - 2))./2) - Dy2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp - Dy1 + Dy2 - 2))./2))./(4.*1.^2.*(cos(pi.*r) + 1))));
        % case 6
        %     % Region F
        %     g = 1./Vin.*sqrt((1/2/pi./r).*  -(Vin.^2.*sin(Dy1.*pi.*r)+Vo.^2.*sin(Dy2.*pi.*r)+2.*Vin.^2.*sin(pi.*r)+2.*Vo.^2.*sin(pi.*r)-2.*Vin.^2.*sin(pi.*r).*cos((Dy1.*pi.*r)/2).^2-2.*Vo.^2.*sin(pi.*r).*cos((Dy2.*pi.*r)/2).^2-2.*Vin.^2.*pi.*r-2.*Vo.^2.*pi.*r+Dy1.*Vin.^2.*pi.*r+Dy2.*Vo.^2.*pi.*r+2.*Vin.^2.*pi.*r.*cos((Dy1.*pi.*r)/2).^2+2.*Vo.^2.*pi.*r.*cos((Dy2.*pi.*r)/2).^2+2.*Vin.^2.*cos(pi.*r).*cos((Dy1.*pi.*r)/2).*sin((Dy1.*pi.*r)/2)+2.*Vo.^2.*cos(pi.*r).*cos((Dy2.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)-2.*Dy1.*Vin.^2.*pi.*r.*cos((Dy1.*pi.*r)/2).^2-2.*Dy2.*Vo.^2.*pi.*r.*cos((Dy2.*pi.*r)/2).^2+Dy1.*Vin.^2.*pi.*r.*cos(pi.*r)+Dy2.*Vo.^2.*pi.*r.*cos(pi.*r)+4.*Vin.*Vo.*sin(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)-2.*Dy1.*Vin.^2.*pi.*r.*cos(pi.*r).*cos((Dy1.*pi.*r)/2).^2-2.*Dy2.*Vo.^2.*pi.*r.*cos(pi.*r).*cos((Dy2.*pi.*r)/2).^2+4.*Vin.*Vo.*cos(pi.*r).*cos(Dp.*pi.*r).*cos((Dy2.*pi.*r)/2).*sin((Dy1.*pi.*r)/2)+4.*Vin.*Vo.*sin(pi.*r).*cos((Dy2.*pi.*r)/2).*sin(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2)+4.*Vin.*Vo.*cos(pi.*r).^2.*cos(Dp.*pi.*r).*cos((Dy2.*pi.*r)/2).*sin((Dy1.*pi.*r)/2)-4.*Vin.*Vo.*cos(pi.*r).^2.*sin(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)+4.*Vin.*Vo.*cos(pi.*r).*sin(pi.*r).*cos(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)+4.*Vin.*Vo.*cos(pi.*r).*sin(pi.*r).*cos((Dy2.*pi.*r)/2).*sin(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2)+4.*Vin.*Vo.*pi.*r.*cos(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)-2.*Dy1.*Vin.^2.*pi.*r.*sin(pi.*r).*cos((Dy1.*pi.*r)/2).*sin((Dy1.*pi.*r)/2)-2.*Dy2.*Vo.^2.*pi.*r.*sin(pi.*r).*cos((Dy2.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)-4.*Dp.*Vin.*Vo.*pi.*r.*cos(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)-2.*Dy1.*Vin.*Vo.*pi.*r.*cos((Dy1.*pi.*r)/2).*sin(Dp.*pi.*r).*sin((Dy2.*pi.*r)/2)-2.*Dy2.*Vin.*Vo.*pi.*r.*cos((Dy2.*pi.*r)/2).*sin(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2)-4.*Vin.*Vo.*pi.*r.*cos(pi.*r).*cos(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)-4.*Vin.*Vo.*pi.*r.*cos(pi.*r).*cos((Dy2.*pi.*r)/2).*sin(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2)+4.*Vin.*Vo.*pi.*r.*sin(pi.*r).*cos(Dp.*pi.*r).*cos((Dy2.*pi.*r)/2).*sin((Dy1.*pi.*r)/2)-4.*Vin.*Vo.*pi.*r.*sin(pi.*r).*sin(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)-4.*Vin.*Vo.*pi.*r.*cos(pi.*r).^2.*cos(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)-4.*Vin.*Vo.*pi.*r.*cos(pi.*r).^2.*cos((Dy2.*pi.*r)/2).*sin(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2)+4.*Dp.*Vin.*Vo.*pi.*r.*cos(pi.*r).^2.*cos(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)+4.*Dp.*Vin.*Vo.*pi.*r.*cos(pi.*r).^2.*cos((Dy2.*pi.*r)/2).*sin(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2)+2.*Dy1.*Vin.*Vo.*pi.*r.*cos(pi.*r).^2.*cos((Dy1.*pi.*r)/2).*sin(Dp.*pi.*r).*sin((Dy2.*pi.*r)/2)+2.*Dy2.*Vin.*Vo.*pi.*r.*cos(pi.*r).^2.*cos(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)+2.*Dy2.*Vin.*Vo.*pi.*r.*cos(pi.*r).^2.*cos((Dy2.*pi.*r)/2).*sin(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2)+4.*Vin.*Vo.*pi.*r.*cos(pi.*r).*sin(pi.*r).*cos(Dp.*pi.*r).*cos((Dy2.*pi.*r)/2).*sin((Dy1.*pi.*r)/2)-4.*Vin.*Vo.*pi.*r.*cos(pi.*r).*sin(pi.*r).*sin(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)-2.*Dy1.*Vin.*Vo.*pi.*r.*cos(pi.*r).*cos(Dp.*pi.*r).*cos((Dy1.*pi.*r)/2).*cos((Dy2.*pi.*r)/2)+4.*Dp.*Vin.*Vo.*pi.*r.*cos(pi.*r).*cos((Dy2.*pi.*r)/2).*sin(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2)-4.*Dp.*Vin.*Vo.*pi.*r.*sin(pi.*r).*cos(Dp.*pi.*r).*cos((Dy2.*pi.*r)/2).*sin((Dy1.*pi.*r)/2)-2.*Dy1.*Vin.*Vo.*pi.*r.*sin(pi.*r).*cos((Dy1.*pi.*r)/2).*cos((Dy2.*pi.*r)/2).*sin(Dp.*pi.*r)+2.*Dy2.*Vin.*Vo.*pi.*r.*cos(pi.*r).*cos(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)+2.*Dy2.*Vin.*Vo.*pi.*r.*sin(pi.*r).*sin(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)-2.*Dy1.*Vin.*Vo.*pi.*r.*cos(pi.*r).^2.*cos(Dp.*pi.*r).*cos((Dy1.*pi.*r)/2).*cos((Dy2.*pi.*r)/2)-4.*Dp.*Vin.*Vo.*pi.*r.*cos(pi.*r).*sin(pi.*r).*cos(Dp.*pi.*r).*cos((Dy2.*pi.*r)/2).*sin((Dy1.*pi.*r)/2)-2.*Dy1.*Vin.*Vo.*pi.*r.*cos(pi.*r).*sin(pi.*r).*cos(Dp.*pi.*r).*cos((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)-2.*Dy1.*Vin.*Vo.*pi.*r.*cos(pi.*r).*sin(pi.*r).*cos((Dy1.*pi.*r)/2).*cos((Dy2.*pi.*r)/2).*sin(Dp.*pi.*r)-2.*Dy2.*Vin.*Vo.*pi.*r.*cos(pi.*r).*sin(pi.*r).*cos(Dp.*pi.*r).*cos((Dy2.*pi.*r)/2).*sin((Dy1.*pi.*r)/2)+4.*Dp.*Vin.*Vo.*pi.*r.*cos(pi.*r).*sin(pi.*r).*sin(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)+2.*Dy2.*Vin.*Vo.*pi.*r.*cos(pi.*r).*sin(pi.*r).*sin(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2))./(2.*1.^2.*(cos(pi.*r)+1)));
        % case 7
        %     % Region G
        %     g = 1./Vin.*sqrt((1/2/pi./r).*  (-(Vin.^2.*sin(Dy1.*pi.*r)+Vo.^2.*sin(Dy2.*pi.*r)+2.*Vin.^2.*sin(pi.*r)+2.*Vo.^2.*sin(pi.*r)-2.*Vin.^2.*sin(pi.*r).*cos((Dy1.*pi.*r)/2).^2-2.*Vo.^2.*sin(pi.*r).*cos((Dy2.*pi.*r)/2).^2-2.*Vin.^2.*pi.*r-2.*Vo.^2.*pi.*r+Dy1.*Vin.^2.*pi.*r+Dy2.*Vo.^2.*pi.*r+2.*Vin.^2.*pi.*r.*cos((Dy1.*pi.*r)/2).^2+2.*Vo.^2.*pi.*r.*cos((Dy2.*pi.*r)/2).^2+2.*Vin.^2.*cos(pi.*r).*cos((Dy1.*pi.*r)/2).*sin((Dy1.*pi.*r)/2)+2.*Vo.^2.*cos(pi.*r).*cos((Dy2.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)-2.*Dy1.*Vin.^2.*pi.*r.*cos((Dy1.*pi.*r)/2).^2-2.*Dy2.*Vo.^2.*pi.*r.*cos((Dy2.*pi.*r)/2).^2+Dy1.*Vin.^2.*pi.*r.*cos(pi.*r)+Dy2.*Vo.^2.*pi.*r.*cos(pi.*r)+4.*Vin.*Vo.*sin(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)-2.*Dy1.*Vin.^2.*pi.*r.*cos(pi.*r).*cos((Dy1.*pi.*r)/2).^2-2.*Dy2.*Vo.^2.*pi.*r.*cos(pi.*r).*cos((Dy2.*pi.*r)/2).^2+4.*Vin.*Vo.*cos(pi.*r).*cos(Dp.*pi.*r).*cos((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)+4.*Vin.*Vo.*sin(pi.*r).*cos((Dy1.*pi.*r)/2).*sin(Dp.*pi.*r).*sin((Dy2.*pi.*r)/2)+4.*Vin.*Vo.*cos(pi.*r).^2.*cos(Dp.*pi.*r).*cos((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)-4.*Vin.*Vo.*cos(pi.*r).^2.*sin(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)+4.*Vin.*Vo.*cos(pi.*r).*sin(pi.*r).*cos(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)+4.*Vin.*Vo.*cos(pi.*r).*sin(pi.*r).*cos((Dy1.*pi.*r)/2).*sin(Dp.*pi.*r).*sin((Dy2.*pi.*r)/2)+4.*Vin.*Vo.*pi.*r.*cos(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)-2.*Dy1.*Vin.^2.*pi.*r.*sin(pi.*r).*cos((Dy1.*pi.*r)/2).*sin((Dy1.*pi.*r)/2)-2.*Dy2.*Vo.^2.*pi.*r.*sin(pi.*r).*cos((Dy2.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)-4.*Dp.*Vin.*Vo.*pi.*r.*cos(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)-2.*Dy1.*Vin.*Vo.*pi.*r.*cos((Dy1.*pi.*r)/2).*sin(Dp.*pi.*r).*sin((Dy2.*pi.*r)/2)-2.*Dy2.*Vin.*Vo.*pi.*r.*cos((Dy2.*pi.*r)/2).*sin(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2)-4.*Vin.*Vo.*pi.*r.*cos(pi.*r).*cos(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)-4.*Vin.*Vo.*pi.*r.*cos(pi.*r).*cos((Dy1.*pi.*r)/2).*sin(Dp.*pi.*r).*sin((Dy2.*pi.*r)/2)+4.*Vin.*Vo.*pi.*r.*sin(pi.*r).*cos(Dp.*pi.*r).*cos((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)-4.*Vin.*Vo.*pi.*r.*sin(pi.*r).*sin(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)-4.*Vin.*Vo.*pi.*r.*cos(pi.*r).^2.*cos(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)-4.*Vin.*Vo.*pi.*r.*cos(pi.*r).^2.*cos((Dy1.*pi.*r)/2).*sin(Dp.*pi.*r).*sin((Dy2.*pi.*r)/2)+4.*Dp.*Vin.*Vo.*pi.*r.*cos(pi.*r).^2.*cos(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)+4.*Dp.*Vin.*Vo.*pi.*r.*cos(pi.*r).^2.*cos((Dy1.*pi.*r)/2).*sin(Dp.*pi.*r).*sin((Dy2.*pi.*r)/2)+2.*Dy1.*Vin.*Vo.*pi.*r.*cos(pi.*r).^2.*cos(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)+2.*Dy1.*Vin.*Vo.*pi.*r.*cos(pi.*r).^2.*cos((Dy1.*pi.*r)/2).*sin(Dp.*pi.*r).*sin((Dy2.*pi.*r)/2)+2.*Dy2.*Vin.*Vo.*pi.*r.*cos(pi.*r).^2.*cos((Dy2.*pi.*r)/2).*sin(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2)+4.*Vin.*Vo.*pi.*r.*cos(pi.*r).*sin(pi.*r).*cos(Dp.*pi.*r).*cos((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)-4.*Vin.*Vo.*pi.*r.*cos(pi.*r).*sin(pi.*r).*sin(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)-2.*Dy2.*Vin.*Vo.*pi.*r.*cos(pi.*r).*cos(Dp.*pi.*r).*cos((Dy1.*pi.*r)/2).*cos((Dy2.*pi.*r)/2)+4.*Dp.*Vin.*Vo.*pi.*r.*cos(pi.*r).*cos((Dy1.*pi.*r)/2).*sin(Dp.*pi.*r).*sin((Dy2.*pi.*r)/2)-4.*Dp.*Vin.*Vo.*pi.*r.*sin(pi.*r).*cos(Dp.*pi.*r).*cos((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)+2.*Dy1.*Vin.*Vo.*pi.*r.*cos(pi.*r).*cos(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)-2.*Dy2.*Vin.*Vo.*pi.*r.*sin(pi.*r).*cos((Dy1.*pi.*r)/2).*cos((Dy2.*pi.*r)/2).*sin(Dp.*pi.*r)+2.*Dy1.*Vin.*Vo.*pi.*r.*sin(pi.*r).*sin(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)-2.*Dy2.*Vin.*Vo.*pi.*r.*cos(pi.*r).^2.*cos(Dp.*pi.*r).*cos((Dy1.*pi.*r)/2).*cos((Dy2.*pi.*r)/2)-4.*Dp.*Vin.*Vo.*pi.*r.*cos(pi.*r).*sin(pi.*r).*cos(Dp.*pi.*r).*cos((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)-2.*Dy1.*Vin.*Vo.*pi.*r.*cos(pi.*r).*sin(pi.*r).*cos(Dp.*pi.*r).*cos((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)-2.*Dy2.*Vin.*Vo.*pi.*r.*cos(pi.*r).*sin(pi.*r).*cos(Dp.*pi.*r).*cos((Dy2.*pi.*r)/2).*sin((Dy1.*pi.*r)/2)-2.*Dy2.*Vin.*Vo.*pi.*r.*cos(pi.*r).*sin(pi.*r).*cos((Dy1.*pi.*r)/2).*cos((Dy2.*pi.*r)/2).*sin(Dp.*pi.*r)+4.*Dp.*Vin.*Vo.*pi.*r.*cos(pi.*r).*sin(pi.*r).*sin(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2)+2.*Dy1.*Vin.*Vo.*pi.*r.*cos(pi.*r).*sin(pi.*r).*sin(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2))./(2.*1.^2.*(cos(pi.*r)+1))));
        % case 8
        %     % Region H
        %     g = 1./Vin.*sqrt((1/2/pi./r).*  (-(2.*Vin.^2.*sin(Dy1.*pi.*r)+2.*Vo.^2.*sin(Dy2.*pi.*r)+2.*Vin.^2.*sin(pi.*r.*(Dy1-1))+2.*Vo.^2.*sin(pi.*r.*(Dy2-1))+2.*Vin.^2.*sin(pi.*r)+2.*Vo.^2.*sin(pi.*r)+2.*Vin.*Vo.*sin((pi.*r.*(Dy1-2.*Dp+Dy2))/2)+2.*Vin.*Vo.*sin((pi.*r.*(Dy1-2.*Dp+Dy2+2))/2)+2.*Vin.*Vo.*sin((pi.*r.*(2.*Dp+Dy1+Dy2-2))/2)+2.*Vin.*Vo.*sin((pi.*r.*(2.*Dp+Dy1+Dy2-4))/2)+2.*Vin.*Vo.*sin((pi.*r.*(2.*Dp+Dy1-Dy2))/2)+2.*Vin.*Vo.*sin((pi.*r.*(2.*Dp-Dy1+Dy2))/2)+2.*Vin.*Vo.*sin((pi.*r.*(2.*Dp+Dy1-Dy2-2))/2)+2.*Vin.*Vo.*sin((pi.*r.*(2.*Dp-Dy1+Dy2-2))/2)-2.*Vin.^2.*pi.*r-2.*Vo.^2.*pi.*r+2.*Vin.^2.*pi.*r.*cos(Dy1.*pi.*r)+2.*Vo.^2.*pi.*r.*cos(Dy2.*pi.*r)-2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(Dy1-2.*Dp+Dy2))/2)+4.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp+Dy1+Dy2-2))/2)+2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp+Dy1+Dy2-4))/2)-2.*Dy1.*Vin.^2.*pi.*r.*cos(Dy1.*pi.*r)-2.*Dy2.*Vo.^2.*pi.*r.*cos(Dy2.*pi.*r)+2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp+Dy1-Dy2))/2)+2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp-Dy1+Dy2))/2)-2.*Dy1.*Vin.^2.*pi.*r.*cos(pi.*r.*(Dy1-1))-2.*Dy2.*Vo.^2.*pi.*r.*cos(pi.*r.*(Dy2-1))+2.*Dp.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(Dy1-2.*Dp+Dy2))/2)+2.*Dp.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(Dy1-2.*Dp+Dy2+2))/2)-2.*Dp.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp+Dy1+Dy2-2))/2)-2.*Dp.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp+Dy1+Dy2-4))/2)-Dy1.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(Dy1-2.*Dp+Dy2))/2)-Dy2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(Dy1-2.*Dp+Dy2))/2)-Dy1.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(Dy1-2.*Dp+Dy2+2))/2)-Dy1.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp+Dy1+Dy2-2))/2)-Dy2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(Dy1-2.*Dp+Dy2+2))/2)-Dy2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp+Dy1+Dy2-2))/2)-Dy1.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp+Dy1+Dy2-4))/2)-Dy2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp+Dy1+Dy2-4))/2)-2.*Dp.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp+Dy1-Dy2))/2)-2.*Dp.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp-Dy1+Dy2))/2)-2.*Dp.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp+Dy1-Dy2-2))/2)-2.*Dp.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp-Dy1+Dy2-2))/2)-Dy1.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp+Dy1-Dy2))/2)+Dy1.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp-Dy1+Dy2))/2)+Dy2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp+Dy1-Dy2))/2)-Dy2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp-Dy1+Dy2))/2)-Dy1.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp+Dy1-Dy2-2))/2)+Dy1.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp-Dy1+Dy2-2))/2)+Dy2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp+Dy1-Dy2-2))/2)-Dy2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp-Dy1+Dy2-2))/2))./(4.*1.^2.*(cos(pi.*r)+1))));
        otherwise
            error('Unknown region');
    end
end


% Define constraint function for ZVS
function [c, ceq] = zvs_function(DoF, vars, Vin, Vo, Z0)
    % Get variables based on DoF
    [Dp, Dy1, Dy2, r] = assign_variables(DoF, vars, Vin, Vo);
    

    region_flag((Dy1>Dy2) & (Dp<(Dy1-Dy2)/2)) = 1;  %  A
    region_flag((Dy2>=Dy1) & (Dp<(Dy2-Dy1)/2)) = 2;  %  B
    region_flag((Dp>=abs(Dy1-Dy2)/2) & (Dp<(Dy2+Dy1)/2) & (Dp<(1-(Dy2+Dy1)/2))) = 3;  %  C
    region_flag((((Dy2+Dy1)/2)<= Dp) & (Dp<(1-(Dy2+Dy1)/2))) = 4;  %  D
    region_flag(((1-(Dy2+Dy1)/2)<= Dp) & (((Dy2+Dy1)/2)>= Dp)) = 5;  %  E
    % region_flag(((1-(Dy2-Dy1)/2)<=Dp) & (Dy2>Dy1)) = 6;  %  F
    % region_flag(((1-(Dy1-Dy2)/2)<=Dp) & (Dy2<=Dy1)) = 7;  %  G
    % region_flag((Dp>=(1-(Dy1+Dy2)/2)) & (Dp>=(Dy1+Dy2)/2) & (Dp<(1-abs(Dy1-Dy2)/2))) = 8;  %  H
    
    % ZVS constraints
    I_zvs_p12_RegionA = -((sin(pi.*r).*(Vo.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2) - Vin.*cos(pi.*r) - Vo.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2) + Vin.*cos(pi.*r.*(Dy1 - 1))))/Z0 + ((cos(pi.*r) + 1).*(Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2) + Vin.*sin(pi.*r) - Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2) + Vin.*sin(pi.*r.*(Dy1 - 1))))/Z0)./(2.*(cos(pi.*r) + 1));
    I_zvs_s12_RegionA = (Vin.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2))/2) + Vo.*sin(pi.*r) + Vin.*sin((r.*pi.*(2.*Dp - Dy1 - Dy2 + 2))/2) + Vo.*sin(pi.*r.*(Dy2 - 1)) + Vin.*cos(pi.*r).*sin((r.*pi.*(2.*Dp - Dy1 - Dy2 + 2))/2) - Vin.*sin(pi.*r).*cos((r.*pi.*(2.*Dp - Dy1 - Dy2 + 2))/2) + Vo.*cos(pi.*r.*(Dy2 - 1)).*sin(pi.*r) + Vo.*sin(pi.*r.*(Dy2 - 1)).*cos(pi.*r) - Vin.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2))/2).*sin(pi.*r) + Vin.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2))/2).*cos(pi.*r))./(2.*Z0.*(cos(pi.*r) + 1));
    I_zvs_p34_RegionA = ((sin(pi.*r).*(2.*Vin - 2.*Vin.*cos((pi.*Dy1.*r)/2).^2 - 2.*Vo.*cos(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vo.*cos((pi.*Dy1.*r)/2).*sin(pi.*Dp.*r).*sin((pi.*Dy2.*r)/2)))/Z0 - ((cos(pi.*r) + 1).*(2.*Vo.*cos(pi.*Dp.*r).*cos((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) - 2.*Vin.*cos((pi.*Dy1.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vo.*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2)))/Z0)./(2.*(cos(pi.*r) + 1));
    I_zvs_s34_RegionA = (((cos(pi.*r) + 1).*(Vin.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2 + 2))/2) - Vin.*sin(pi.*r) - Vo.*sin(pi.*Dy2.*r) + Vin.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2))/2)))/Z0 + (sin(pi.*r).*(Vin - Vo - Vin.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2 + 2))/2) + Vin.*cos(pi.*r) + Vo.*cos(pi.*Dy2.*r) - Vin.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));

    I_zvs_p12_RegionB = (Vin.*sin(pi.*r).*cos(pi.*Dy1.*r) - Vin.*sin(pi.*Dy1.*r) - Vin.*cos(pi.*r).*sin(pi.*Dy1.*r) - Vin.*sin(pi.*r) + 2.*Vo.*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vo.*cos((pi.*Dy1.*r)/2).*cos((pi.*Dy2.*r)/2).*sin(pi.*Dp.*r) + 2.*Vo.*cos(pi.*r).*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vo.*cos(pi.*r).*cos((pi.*Dy1.*r)/2).*cos((pi.*Dy2.*r)/2).*sin(pi.*Dp.*r) + 2.*Vo.*sin(pi.*r).*cos(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vo.*sin(pi.*r).*cos((pi.*Dy1.*r)/2).*sin(pi.*Dp.*r).*sin((pi.*Dy2.*r)/2))./(2.*Z0.*(cos(pi.*r) + 1));
    I_zvs_s12_RegionB = -(Vo.*sin(pi.*r).*cos((pi.*Dy2.*r)/2).^2 - (Vo.*sin(pi.*Dy2.*r))/2 - Vo.*sin(pi.*r) + Vin.*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + Vin.*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) - Vo.*cos(pi.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy2.*r)/2) + Vin.*cos(pi.*r).*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + Vin.*cos(pi.*r).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + Vin.*sin(pi.*r).*cos(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) - Vin.*sin(pi.*r).*cos((pi.*Dy2.*r)/2).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2))./(Z0.*(cos(pi.*r) + 1));
    I_zvs_p34_RegionB = (Vin.*sin(pi.*r) + Vin.*sin(pi.*Dy1.*r) + Vin.*cos(pi.*r).*sin(pi.*Dy1.*r) - Vin.*sin(pi.*r).*cos(pi.*Dy1.*r) - 2.*Vo.*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vo.*cos((pi.*Dy1.*r)/2).*cos((pi.*Dy2.*r)/2).*sin(pi.*Dp.*r) - 2.*Vo.*cos(pi.*r).*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vo.*cos(pi.*r).*cos((pi.*Dy1.*r)/2).*cos((pi.*Dy2.*r)/2).*sin(pi.*Dp.*r) - 2.*Vo.*sin(pi.*r).*cos(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vo.*sin(pi.*r).*cos((pi.*Dy1.*r)/2).*sin(pi.*Dp.*r).*sin((pi.*Dy2.*r)/2))./(2.*Z0.*(cos(pi.*r) + 1));
    I_zvs_s34_RegionB = -(Vo.*sin(pi.*r) + (Vo.*sin(pi.*Dy2.*r))/2 - Vo.*sin(pi.*r).*cos((pi.*Dy2.*r)/2).^2 - Vin.*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + Vin.*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + Vo.*cos(pi.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy2.*r)/2) - Vin.*cos(pi.*r).*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + Vin.*cos(pi.*r).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) - Vin.*sin(pi.*r).*cos(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) - Vin.*sin(pi.*r).*cos((pi.*Dy2.*r)/2).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2))./(Z0.*(cos(pi.*r) + 1));

    I_zvs_p12_RegionC = -((sin(pi.*r).*(Vo.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2) - Vin.*cos(pi.*r) - Vo.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2) + Vin.*cos(pi.*r.*(Dy1 - 1))))/Z0 + ((cos(pi.*r) + 1).*(Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2) + Vin.*sin(pi.*r) - Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2) + Vin.*sin(pi.*r.*(Dy1 - 1))))/Z0)./(2.*(cos(pi.*r) + 1));
    I_zvs_s12_RegionC = (Vin.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2))/2) + Vo.*sin(pi.*r) + Vin.*sin((r.*pi.*(2.*Dp - Dy1 - Dy2 + 2))/2) + Vo.*sin(pi.*r.*(Dy2 - 1)) + Vin.*cos(pi.*r).*sin((r.*pi.*(2.*Dp - Dy1 - Dy2 + 2))/2) - Vin.*sin(pi.*r).*cos((r.*pi.*(2.*Dp - Dy1 - Dy2 + 2))/2) + Vo.*cos(pi.*r.*(Dy2 - 1)).*sin(pi.*r) + Vo.*sin(pi.*r.*(Dy2 - 1)).*cos(pi.*r) - Vin.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2))/2).*sin(pi.*r) + Vin.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2))/2).*cos(pi.*r))./(2.*Z0.*(cos(pi.*r) + 1));
    I_zvs_p34_RegionC = (((cos(pi.*r) + 1).*(Vo.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) + Vo.*sin(pi.*r) + Vin.*sin(pi.*Dy1.*r) - Vo.*sin((r.*pi.*(Dy1 - 2.*Dp + Dy2))/2)))/Z0 + (sin(pi.*r).*(Vin - Vo + Vo.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) - Vo.*cos(pi.*r) - Vin.*cos(pi.*Dy1.*r) + Vo.*cos((r.*pi.*(Dy1 - 2.*Dp + Dy2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));
    I_zvs_s34_RegionC = ((sin(pi.*r).*(2.*Vo.*cos((pi.*Dy2.*r)/2).^2 - 2.*Vo + 2.*Vin.*cos(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vin.*cos((pi.*Dy2.*r)/2).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2)))/Z0 - ((cos(pi.*r) + 1).*(2.*Vo.*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy2.*r)/2) - 2.*Vin.*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vin.*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2)))/Z0)./(2.*(cos(pi.*r) + 1));

    I_zvs_p12_RegionD = -((sin(pi.*r).*(Vo.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2) - Vin.*cos(pi.*r) - Vo.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2) + Vin.*cos(pi.*r.*(Dy1 - 1))))/Z0 + ((cos(pi.*r) + 1).*(Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2) + Vin.*sin(pi.*r) - Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2) + Vin.*sin(pi.*r.*(Dy1 - 1))))/Z0)./(2.*(cos(pi.*r) + 1));
    I_zvs_s12_RegionD = (Vo.*sin(pi.*r) + Vo.*sin(pi.*Dy2.*r) + Vo.*cos(pi.*r).*sin(pi.*Dy2.*r) - Vo.*sin(pi.*r).*cos(pi.*Dy2.*r) + 2.*Vin.*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vin.*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vin.*cos(pi.*r).*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vin.*cos(pi.*r).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) - 2.*Vin.*sin(pi.*r).*cos(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vin.*sin(pi.*r).*cos((pi.*Dy2.*r)/2).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2))./(2.*Z0.*(cos(pi.*r) + 1));
    I_zvs_p34_RegionD = ((sin(pi.*r).*(Vin + Vo.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) - Vin.*cos(pi.*Dy1.*r) - Vo.*cos((r.*pi.*(Dy1 - 2.*Dp + Dy2 + 2))/2)))/Z0 + ((cos(pi.*r) + 1).*(Vo.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) + Vin.*sin(pi.*Dy1.*r) + Vo.*sin((r.*pi.*(Dy1 - 2.*Dp + Dy2 + 2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));
    I_zvs_s34_RegionD = ((sin(pi.*r).*(2.*Vo.*cos((pi.*Dy2.*r)/2).^2 - 2.*Vo + 2.*Vin.*cos(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vin.*cos((pi.*Dy2.*r)/2).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2)))/Z0 - ((cos(pi.*r) + 1).*(2.*Vo.*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy2.*r)/2) - 2.*Vin.*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vin.*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2)))/Z0)./(2.*(cos(pi.*r) + 1));

    I_zvs_p12_RegionE = -(Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2) + Vin.*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2) + Vin.*sin(pi.*r.*(Dy1 - 1)) + Vo.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2).*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2).*cos(pi.*r) + Vin.*cos(pi.*r.*(Dy1 - 1)).*sin(pi.*r) + Vin.*sin(pi.*r.*(Dy1 - 1)).*cos(pi.*r) + Vo.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2).*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2).*cos(pi.*r))./(2.*Z0.*(cos(pi.*r) + 1));
    I_zvs_s12_RegionE = (Vin.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2))/2) + Vo.*sin(pi.*r) + Vin.*sin((r.*pi.*(2.*Dp - Dy1 - Dy2 + 2))/2) + Vo.*sin(pi.*r.*(Dy2 - 1)) + Vin.*cos(pi.*r).*sin((r.*pi.*(2.*Dp - Dy1 - Dy2 + 2))/2) - Vin.*sin(pi.*r).*cos((r.*pi.*(2.*Dp - Dy1 - Dy2 + 2))/2) + Vo.*cos(pi.*r.*(Dy2 - 1)).*sin(pi.*r) + Vo.*sin(pi.*r.*(Dy2 - 1)).*cos(pi.*r) - Vin.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2))/2).*sin(pi.*r) + Vin.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2))/2).*cos(pi.*r))./(2.*Z0.*(cos(pi.*r) + 1));
    I_zvs_p34_RegionE =  (((cos(pi.*r) + 1).*(Vo.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) + Vo.*sin(pi.*r) + Vin.*sin(pi.*Dy1.*r) - Vo.*sin((r.*pi.*(Dy1 - 2.*Dp + Dy2))/2)))/Z0 + (sin(pi.*r).*(Vin - Vo + Vo.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) - Vo.*cos(pi.*r) - Vin.*cos(pi.*Dy1.*r) + Vo.*cos((r.*pi.*(Dy1 - 2.*Dp + Dy2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));
    I_zvs_s34_RegionE = -(((cos(pi.*r) + 1).*(Vin.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2))/2) - Vin.*sin(pi.*r) + Vo.*sin(pi.*Dy2.*r) + Vin.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2)))/Z0 + (sin(pi.*r).*(Vin + Vo - Vin.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2))/2) + Vin.*cos(pi.*r) - Vo.*cos(pi.*Dy2.*r) - Vin.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));

    % I_zvs_p12_RegionF = -(Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2) + Vin.*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2) + Vin.*sin(pi.*r.*(Dy1 - 1)) + Vo.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2).*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2).*cos(pi.*r) + Vin.*cos(pi.*r.*(Dy1 - 1)).*sin(pi.*r) + Vin.*sin(pi.*r.*(Dy1 - 1)).*cos(pi.*r) + Vo.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2).*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2).*cos(pi.*r))./(2.*Z0.*(cos(pi.*r) + 1));
    % I_zvs_s12_RegionF = (Vo.*sin(pi.*r) + Vo.*sin(pi.*Dy2.*r) + Vo.*cos(pi.*r).*sin(pi.*Dy2.*r) - Vo.*sin(pi.*r).*cos(pi.*Dy2.*r) + 2.*Vin.*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vin.*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vin.*cos(pi.*r).*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vin.*cos(pi.*r).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) - 2.*Vin.*sin(pi.*r).*cos(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vin.*sin(pi.*r).*cos((pi.*Dy2.*r)/2).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2))./(2.*Z0.*(cos(pi.*r) + 1));
    % I_zvs_p34_RegionF = -(((cos(pi.*r) + 1).*(Vo.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2 - 4))/2) + Vo.*sin(pi.*r) - Vin.*sin(pi.*Dy1.*r) - Vo.*sin((r.*pi.*(Dy1 - 2.*Dp + Dy2 + 2))/2)))/Z0 - (sin(pi.*r).*(Vin + Vo - Vo.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2 - 4))/2) + Vo.*cos(pi.*r) - Vin.*cos(pi.*Dy1.*r) - Vo.*cos((r.*pi.*(Dy1 - 2.*Dp + Dy2 + 2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));
    % I_zvs_s34_RegionF = -((sin(pi.*r).*(Vo + Vin.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) - Vo.*cos(pi.*Dy2.*r) - Vin.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2)))/Z0 + ((cos(pi.*r) + 1).*(Vo.*sin(pi.*Dy2.*r) - Vin.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) + Vin.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));
    % 
    % I_zvs_p12_RegionG = ((sin(pi*r)*(Vo*cos((r*pi*(2*Dp + Dy1 - Dy2 - 4))/2) + Vin*cos(pi*r) - Vo*cos((r*pi*(2*Dp + Dy1 + Dy2 - 4))/2) - Vin*cos(pi*r*(Dy1 - 1))))/Z0 - ((cos(pi*r) + 1)*(Vin*sin(pi*r) - Vo*sin((r*pi*(2*Dp + Dy1 - Dy2 - 4))/2) + Vo*sin((r*pi*(2*Dp + Dy1 + Dy2 - 4))/2) + Vin*sin(pi*r*(Dy1 - 1))))/Z0)/(2*(cos(pi*r) + 1));
    % I_zvs_s12_RegionG = (Vo*sin(pi*r) - Vin*sin((r*pi*(2*Dp + Dy1 - Dy2 - 2))/2) + Vin*sin((r*pi*(Dy1 - 2*Dp + Dy2))/2) + Vo*sin(pi*r*(Dy2 - 1)) + Vin*cos((r*pi*(Dy1 - 2*Dp + Dy2))/2)*sin(pi*r) + Vin*sin((r*pi*(Dy1 - 2*Dp + Dy2))/2)*cos(pi*r) + Vo*cos(pi*r*(Dy2 - 1))*sin(pi*r) + Vo*sin(pi*r*(Dy2 - 1))*cos(pi*r) + Vin*cos((r*pi*(2*Dp + Dy1 - Dy2 - 2))/2)*sin(pi*r) - Vin*sin((r*pi*(2*Dp + Dy1 - Dy2 - 2))/2)*cos(pi*r))/(2*Z0*(cos(pi*r) + 1));
    % I_zvs_p34_RegionG = ((sin(pi*r)*(Vin + Vo*cos((r*pi*(2*Dp - Dy1 + Dy2 - 2))/2) - Vin*cos(pi*Dy1*r) - Vo*cos((r*pi*(Dy1 - 2*Dp + Dy2 + 2))/2)))/Z0 + ((cos(pi*r) + 1)*(Vo*sin((r*pi*(2*Dp - Dy1 + Dy2 - 2))/2) + Vin*sin(pi*Dy1*r) + Vo*sin((r*pi*(Dy1 - 2*Dp + Dy2 + 2))/2)))/Z0)/(2*(cos(pi*r) + 1));
    % I_zvs_s34_RegionG = -(((cos(pi*r) + 1)*(Vin*sin((r*pi*(2*Dp - Dy1 + Dy2))/2) - Vin*sin(pi*r) + Vo*sin(pi*Dy2*r) + Vin*sin((r*pi*(2*Dp + Dy1 + Dy2 - 2))/2)))/Z0 + (sin(pi*r)*(Vin + Vo - Vin*cos((r*pi*(2*Dp - Dy1 + Dy2))/2) + Vin*cos(pi*r) - Vo*cos(pi*Dy2*r) - Vin*cos((r*pi*(2*Dp + Dy1 + Dy2 - 2))/2)))/Z0)/(2*(cos(pi*r) + 1));
    % 
    % I_zvs_p12_RegionH = -(Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2) + Vin.*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2) + Vin.*sin(pi.*r.*(Dy1 - 1)) + Vo.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2).*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2).*cos(pi.*r) + Vin.*cos(pi.*r.*(Dy1 - 1)).*sin(pi.*r) + Vin.*sin(pi.*r.*(Dy1 - 1)).*cos(pi.*r) + Vo.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2).*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2).*cos(pi.*r))./(2.*Z0.*(cos(pi.*r) + 1));
    % I_zvs_s12_RegionH = (Vo.*sin(pi.*r) + Vo.*sin(pi.*Dy2.*r) + Vo.*cos(pi.*r).*sin(pi.*Dy2.*r) - Vo.*sin(pi.*r).*cos(pi.*Dy2.*r) + 2.*Vin.*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vin.*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vin.*cos(pi.*r).*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vin.*cos(pi.*r).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) - 2.*Vin.*sin(pi.*r).*cos(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vin.*sin(pi.*r).*cos((pi.*Dy2.*r)/2).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2))./(2.*Z0.*(cos(pi.*r) + 1));
    % I_zvs_p34_RegionH = ((sin(pi.*r).*(Vin + Vo.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) - Vin.*cos(pi.*Dy1.*r) - Vo.*cos((r.*pi.*(Dy1 - 2.*Dp + Dy2 + 2))/2)))/Z0 + ((cos(pi.*r) + 1).*(Vo.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) + Vin.*sin(pi.*Dy1.*r) + Vo.*sin((r.*pi.*(Dy1 - 2.*Dp + Dy2 + 2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));
    % I_zvs_s34_RegionH = -(((cos(pi.*r) + 1).*(Vin.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2))/2) - Vin.*sin(pi.*r) + Vo.*sin(pi.*Dy2.*r) + Vin.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2)))/Z0 + (sin(pi.*r).*(Vin + Vo - Vin.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2))/2) + Vin.*cos(pi.*r) - Vo.*cos(pi.*Dy2.*r) - Vin.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));

    % % Define inequality constraints (c <= 0)
    % c = [ I_zvs_p12_RegionE;    % I_zvs_p12_RegionE <= 0
    %      -I_zvs_s12_RegionE;    % -I_zvs_s12_RegionE <= 0 (equivalent to I_zvs_s12_RegionE >= 0)
    %      -I_zvs_p34_RegionE;    % -I_zvs_p34_RegionE <= 0 (equivalent to I_zvs_p34_RegionE >= 0)
    %      I_zvs_s34_RegionE];    % I_zvs_s34_RegionE <= 0

    switch region_flag
        case 1
            % Region A
            c = [I_zvs_p12_RegionA; -I_zvs_s12_RegionA; -I_zvs_p34_RegionA; I_zvs_s34_RegionA];
        case 2
            % Region B
            c = [I_zvs_p12_RegionB; -I_zvs_s12_RegionB; -I_zvs_p34_RegionB; I_zvs_s34_RegionB];
        case 3
            % Region C
            c = [I_zvs_p12_RegionC; -I_zvs_s12_RegionC; -I_zvs_p34_RegionC; I_zvs_s34_RegionC];
        case 4
            % Region D
            c = [I_zvs_p12_RegionD; -I_zvs_s12_RegionD; -I_zvs_p34_RegionD; I_zvs_s34_RegionD];
        case 5
            % Region E
            c = [I_zvs_p12_RegionE; -I_zvs_s12_RegionE; -I_zvs_p34_RegionE; I_zvs_s34_RegionE];
        % case 6
        %     % Region F
        %     c = [I_zvs_p12_RegionF; -I_zvs_s12_RegionF; -I_zvs_p34_RegionF; I_zvs_s34_RegionF];
        % case 7
        %     % Region G
        %     c = [I_zvs_p12_RegionG; -I_zvs_s12_RegionG; -I_zvs_p34_RegionG; I_zvs_s34_RegionG];
        % case 8
        %     % Region H
        %     c = [I_zvs_p12_RegionH; -I_zvs_s12_RegionH; -I_zvs_p34_RegionH; I_zvs_s34_RegionH];
        otherwise
            error('Unknown region');
    end

    % Define equality constraint (ceq = 0)
    ceq = [];     % PoN should equal PoN_target
end

% Define Ioff function
function [Ip12, Is12, Ip34, Is34, ILmax] = h_function(DoF, vars, Vin, Vo)
    Z0 = 1;
    % Get variables based on DoF
    [Dp, Dy1, Dy2, r] = assign_variables(DoF, vars, Vin, Vo);
    

    region_flag((Dy1>Dy2) & (Dp<(Dy1-Dy2)/2)) = 1;  %  A
    region_flag((Dy2>=Dy1) & (Dp<(Dy2-Dy1)/2)) = 2;  %  B
    region_flag((Dp>=abs(Dy1-Dy2)/2) & (Dp<(Dy2+Dy1)/2) & (Dp<(1-(Dy2+Dy1)/2))) = 3;  %  C
    region_flag((((Dy2+Dy1)/2)<= Dp) & (Dp<(1-(Dy2+Dy1)/2))) = 4;  %  D
    region_flag(((1-(Dy2+Dy1)/2)<= Dp) & (((Dy2+Dy1)/2)>= Dp)) = 5;  %  E
    % region_flag(((1-(Dy2-Dy1)/2)<=Dp) & (Dy2>Dy1)) = 6;  %  F
    % region_flag(((1-(Dy1-Dy2)/2)<=Dp) & (Dy2<=Dy1)) = 7;  %  G
    % region_flag((Dp>=(1-(Dy1+Dy2)/2)) & (Dp>=(Dy1+Dy2)/2) & (Dp<(1-abs(Dy1-Dy2)/2))) = 8;  %  H

        % ZVS constraints
        I_zvs_p12_RegionA = -((sin(pi.*r).*(Vo.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2) - Vin.*cos(pi.*r) - Vo.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2) + Vin.*cos(pi.*r.*(Dy1 - 1))))/Z0 + ((cos(pi.*r) + 1).*(Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2) + Vin.*sin(pi.*r) - Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2) + Vin.*sin(pi.*r.*(Dy1 - 1))))/Z0)./(2.*(cos(pi.*r) + 1));
        I_zvs_s12_RegionA = (Vin.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2))/2) + Vo.*sin(pi.*r) + Vin.*sin((r.*pi.*(2.*Dp - Dy1 - Dy2 + 2))/2) + Vo.*sin(pi.*r.*(Dy2 - 1)) + Vin.*cos(pi.*r).*sin((r.*pi.*(2.*Dp - Dy1 - Dy2 + 2))/2) - Vin.*sin(pi.*r).*cos((r.*pi.*(2.*Dp - Dy1 - Dy2 + 2))/2) + Vo.*cos(pi.*r.*(Dy2 - 1)).*sin(pi.*r) + Vo.*sin(pi.*r.*(Dy2 - 1)).*cos(pi.*r) - Vin.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2))/2).*sin(pi.*r) + Vin.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2))/2).*cos(pi.*r))./(2.*Z0.*(cos(pi.*r) + 1));
        I_zvs_p34_RegionA = ((sin(pi.*r).*(2.*Vin - 2.*Vin.*cos((pi.*Dy1.*r)/2).^2 - 2.*Vo.*cos(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vo.*cos((pi.*Dy1.*r)/2).*sin(pi.*Dp.*r).*sin((pi.*Dy2.*r)/2)))/Z0 - ((cos(pi.*r) + 1).*(2.*Vo.*cos(pi.*Dp.*r).*cos((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) - 2.*Vin.*cos((pi.*Dy1.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vo.*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2)))/Z0)./(2.*(cos(pi.*r) + 1));
        I_zvs_s34_RegionA = (((cos(pi.*r) + 1).*(Vin.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2 + 2))/2) - Vin.*sin(pi.*r) - Vo.*sin(pi.*Dy2.*r) + Vin.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2))/2)))/Z0 + (sin(pi.*r).*(Vin - Vo - Vin.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2 + 2))/2) + Vin.*cos(pi.*r) + Vo.*cos(pi.*Dy2.*r) - Vin.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));
    
        I_zvs_p12_RegionB = (Vin.*sin(pi.*r).*cos(pi.*Dy1.*r) - Vin.*sin(pi.*Dy1.*r) - Vin.*cos(pi.*r).*sin(pi.*Dy1.*r) - Vin.*sin(pi.*r) + 2.*Vo.*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vo.*cos((pi.*Dy1.*r)/2).*cos((pi.*Dy2.*r)/2).*sin(pi.*Dp.*r) + 2.*Vo.*cos(pi.*r).*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vo.*cos(pi.*r).*cos((pi.*Dy1.*r)/2).*cos((pi.*Dy2.*r)/2).*sin(pi.*Dp.*r) + 2.*Vo.*sin(pi.*r).*cos(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vo.*sin(pi.*r).*cos((pi.*Dy1.*r)/2).*sin(pi.*Dp.*r).*sin((pi.*Dy2.*r)/2))./(2.*Z0.*(cos(pi.*r) + 1));
        I_zvs_s12_RegionB = -(Vo.*sin(pi.*r).*cos((pi.*Dy2.*r)/2).^2 - (Vo.*sin(pi.*Dy2.*r))/2 - Vo.*sin(pi.*r) + Vin.*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + Vin.*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) - Vo.*cos(pi.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy2.*r)/2) + Vin.*cos(pi.*r).*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + Vin.*cos(pi.*r).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + Vin.*sin(pi.*r).*cos(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) - Vin.*sin(pi.*r).*cos((pi.*Dy2.*r)/2).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2))./(Z0.*(cos(pi.*r) + 1));
        I_zvs_p34_RegionB = (Vin.*sin(pi.*r) + Vin.*sin(pi.*Dy1.*r) + Vin.*cos(pi.*r).*sin(pi.*Dy1.*r) - Vin.*sin(pi.*r).*cos(pi.*Dy1.*r) - 2.*Vo.*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vo.*cos((pi.*Dy1.*r)/2).*cos((pi.*Dy2.*r)/2).*sin(pi.*Dp.*r) - 2.*Vo.*cos(pi.*r).*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vo.*cos(pi.*r).*cos((pi.*Dy1.*r)/2).*cos((pi.*Dy2.*r)/2).*sin(pi.*Dp.*r) - 2.*Vo.*sin(pi.*r).*cos(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vo.*sin(pi.*r).*cos((pi.*Dy1.*r)/2).*sin(pi.*Dp.*r).*sin((pi.*Dy2.*r)/2))./(2.*Z0.*(cos(pi.*r) + 1));
        I_zvs_s34_RegionB = -(Vo.*sin(pi.*r) + (Vo.*sin(pi.*Dy2.*r))/2 - Vo.*sin(pi.*r).*cos((pi.*Dy2.*r)/2).^2 - Vin.*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + Vin.*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + Vo.*cos(pi.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy2.*r)/2) - Vin.*cos(pi.*r).*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + Vin.*cos(pi.*r).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) - Vin.*sin(pi.*r).*cos(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) - Vin.*sin(pi.*r).*cos((pi.*Dy2.*r)/2).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2))./(Z0.*(cos(pi.*r) + 1));
    
        I_zvs_p12_RegionC = -((sin(pi.*r).*(Vo.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2) - Vin.*cos(pi.*r) - Vo.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2) + Vin.*cos(pi.*r.*(Dy1 - 1))))/Z0 + ((cos(pi.*r) + 1).*(Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2) + Vin.*sin(pi.*r) - Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2) + Vin.*sin(pi.*r.*(Dy1 - 1))))/Z0)./(2.*(cos(pi.*r) + 1));
        I_zvs_s12_RegionC = (Vin.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2))/2) + Vo.*sin(pi.*r) + Vin.*sin((r.*pi.*(2.*Dp - Dy1 - Dy2 + 2))/2) + Vo.*sin(pi.*r.*(Dy2 - 1)) + Vin.*cos(pi.*r).*sin((r.*pi.*(2.*Dp - Dy1 - Dy2 + 2))/2) - Vin.*sin(pi.*r).*cos((r.*pi.*(2.*Dp - Dy1 - Dy2 + 2))/2) + Vo.*cos(pi.*r.*(Dy2 - 1)).*sin(pi.*r) + Vo.*sin(pi.*r.*(Dy2 - 1)).*cos(pi.*r) - Vin.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2))/2).*sin(pi.*r) + Vin.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2))/2).*cos(pi.*r))./(2.*Z0.*(cos(pi.*r) + 1));
        I_zvs_p34_RegionC = (((cos(pi.*r) + 1).*(Vo.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) + Vo.*sin(pi.*r) + Vin.*sin(pi.*Dy1.*r) - Vo.*sin((r.*pi.*(Dy1 - 2.*Dp + Dy2))/2)))/Z0 + (sin(pi.*r).*(Vin - Vo + Vo.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) - Vo.*cos(pi.*r) - Vin.*cos(pi.*Dy1.*r) + Vo.*cos((r.*pi.*(Dy1 - 2.*Dp + Dy2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));
        I_zvs_s34_RegionC = ((sin(pi.*r).*(2.*Vo.*cos((pi.*Dy2.*r)/2).^2 - 2.*Vo + 2.*Vin.*cos(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vin.*cos((pi.*Dy2.*r)/2).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2)))/Z0 - ((cos(pi.*r) + 1).*(2.*Vo.*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy2.*r)/2) - 2.*Vin.*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vin.*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2)))/Z0)./(2.*(cos(pi.*r) + 1));
    
        I_zvs_p12_RegionD = -((sin(pi.*r).*(Vo.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2) - Vin.*cos(pi.*r) - Vo.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2) + Vin.*cos(pi.*r.*(Dy1 - 1))))/Z0 + ((cos(pi.*r) + 1).*(Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2) + Vin.*sin(pi.*r) - Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2) + Vin.*sin(pi.*r.*(Dy1 - 1))))/Z0)./(2.*(cos(pi.*r) + 1));
        I_zvs_s12_RegionD = (Vo.*sin(pi.*r) + Vo.*sin(pi.*Dy2.*r) + Vo.*cos(pi.*r).*sin(pi.*Dy2.*r) - Vo.*sin(pi.*r).*cos(pi.*Dy2.*r) + 2.*Vin.*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vin.*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vin.*cos(pi.*r).*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vin.*cos(pi.*r).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) - 2.*Vin.*sin(pi.*r).*cos(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vin.*sin(pi.*r).*cos((pi.*Dy2.*r)/2).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2))./(2.*Z0.*(cos(pi.*r) + 1));
        I_zvs_p34_RegionD = ((sin(pi.*r).*(Vin + Vo.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) - Vin.*cos(pi.*Dy1.*r) - Vo.*cos((r.*pi.*(Dy1 - 2.*Dp + Dy2 + 2))/2)))/Z0 + ((cos(pi.*r) + 1).*(Vo.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) + Vin.*sin(pi.*Dy1.*r) + Vo.*sin((r.*pi.*(Dy1 - 2.*Dp + Dy2 + 2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));
        I_zvs_s34_RegionD = ((sin(pi.*r).*(2.*Vo.*cos((pi.*Dy2.*r)/2).^2 - 2.*Vo + 2.*Vin.*cos(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vin.*cos((pi.*Dy2.*r)/2).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2)))/Z0 - ((cos(pi.*r) + 1).*(2.*Vo.*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy2.*r)/2) - 2.*Vin.*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vin.*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2)))/Z0)./(2.*(cos(pi.*r) + 1));
    
        I_zvs_p12_RegionE = -(Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2) + Vin.*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2) + Vin.*sin(pi.*r.*(Dy1 - 1)) + Vo.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2).*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2).*cos(pi.*r) + Vin.*cos(pi.*r.*(Dy1 - 1)).*sin(pi.*r) + Vin.*sin(pi.*r.*(Dy1 - 1)).*cos(pi.*r) + Vo.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2).*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2).*cos(pi.*r))./(2.*Z0.*(cos(pi.*r) + 1));
        I_zvs_s12_RegionE = (Vin.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2))/2) + Vo.*sin(pi.*r) + Vin.*sin((r.*pi.*(2.*Dp - Dy1 - Dy2 + 2))/2) + Vo.*sin(pi.*r.*(Dy2 - 1)) + Vin.*cos(pi.*r).*sin((r.*pi.*(2.*Dp - Dy1 - Dy2 + 2))/2) - Vin.*sin(pi.*r).*cos((r.*pi.*(2.*Dp - Dy1 - Dy2 + 2))/2) + Vo.*cos(pi.*r.*(Dy2 - 1)).*sin(pi.*r) + Vo.*sin(pi.*r.*(Dy2 - 1)).*cos(pi.*r) - Vin.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2))/2).*sin(pi.*r) + Vin.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2))/2).*cos(pi.*r))./(2.*Z0.*(cos(pi.*r) + 1));
        I_zvs_p34_RegionE =  (((cos(pi.*r) + 1).*(Vo.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) + Vo.*sin(pi.*r) + Vin.*sin(pi.*Dy1.*r) - Vo.*sin((r.*pi.*(Dy1 - 2.*Dp + Dy2))/2)))/Z0 + (sin(pi.*r).*(Vin - Vo + Vo.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) - Vo.*cos(pi.*r) - Vin.*cos(pi.*Dy1.*r) + Vo.*cos((r.*pi.*(Dy1 - 2.*Dp + Dy2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));
        I_zvs_s34_RegionE = -(((cos(pi.*r) + 1).*(Vin.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2))/2) - Vin.*sin(pi.*r) + Vo.*sin(pi.*Dy2.*r) + Vin.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2)))/Z0 + (sin(pi.*r).*(Vin + Vo - Vin.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2))/2) + Vin.*cos(pi.*r) - Vo.*cos(pi.*Dy2.*r) - Vin.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));
    
        % I_zvs_p12_RegionF = -(Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2) + Vin.*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2) + Vin.*sin(pi.*r.*(Dy1 - 1)) + Vo.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2).*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2).*cos(pi.*r) + Vin.*cos(pi.*r.*(Dy1 - 1)).*sin(pi.*r) + Vin.*sin(pi.*r.*(Dy1 - 1)).*cos(pi.*r) + Vo.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2).*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2).*cos(pi.*r))./(2.*Z0.*(cos(pi.*r) + 1));
        % I_zvs_s12_RegionF = (Vo.*sin(pi.*r) + Vo.*sin(pi.*Dy2.*r) + Vo.*cos(pi.*r).*sin(pi.*Dy2.*r) - Vo.*sin(pi.*r).*cos(pi.*Dy2.*r) + 2.*Vin.*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vin.*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vin.*cos(pi.*r).*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vin.*cos(pi.*r).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) - 2.*Vin.*sin(pi.*r).*cos(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vin.*sin(pi.*r).*cos((pi.*Dy2.*r)/2).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2))./(2.*Z0.*(cos(pi.*r) + 1));
        % I_zvs_p34_RegionF = -(((cos(pi.*r) + 1).*(Vo.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2 - 4))/2) + Vo.*sin(pi.*r) - Vin.*sin(pi.*Dy1.*r) - Vo.*sin((r.*pi.*(Dy1 - 2.*Dp + Dy2 + 2))/2)))/Z0 - (sin(pi.*r).*(Vin + Vo - Vo.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2 - 4))/2) + Vo.*cos(pi.*r) - Vin.*cos(pi.*Dy1.*r) - Vo.*cos((r.*pi.*(Dy1 - 2.*Dp + Dy2 + 2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));
        % I_zvs_s34_RegionF = -((sin(pi.*r).*(Vo + Vin.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) - Vo.*cos(pi.*Dy2.*r) - Vin.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2)))/Z0 + ((cos(pi.*r) + 1).*(Vo.*sin(pi.*Dy2.*r) - Vin.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) + Vin.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));
        % 
        % I_zvs_p12_RegionG = ((sin(pi*r)*(Vo*cos((r*pi*(2*Dp + Dy1 - Dy2 - 4))/2) + Vin*cos(pi*r) - Vo*cos((r*pi*(2*Dp + Dy1 + Dy2 - 4))/2) - Vin*cos(pi*r*(Dy1 - 1))))/Z0 - ((cos(pi*r) + 1)*(Vin*sin(pi*r) - Vo*sin((r*pi*(2*Dp + Dy1 - Dy2 - 4))/2) + Vo*sin((r*pi*(2*Dp + Dy1 + Dy2 - 4))/2) + Vin*sin(pi*r*(Dy1 - 1))))/Z0)/(2*(cos(pi*r) + 1));
        % I_zvs_s12_RegionG = (Vo*sin(pi*r) - Vin*sin((r*pi*(2*Dp + Dy1 - Dy2 - 2))/2) + Vin*sin((r*pi*(Dy1 - 2*Dp + Dy2))/2) + Vo*sin(pi*r*(Dy2 - 1)) + Vin*cos((r*pi*(Dy1 - 2*Dp + Dy2))/2)*sin(pi*r) + Vin*sin((r*pi*(Dy1 - 2*Dp + Dy2))/2)*cos(pi*r) + Vo*cos(pi*r*(Dy2 - 1))*sin(pi*r) + Vo*sin(pi*r*(Dy2 - 1))*cos(pi*r) + Vin*cos((r*pi*(2*Dp + Dy1 - Dy2 - 2))/2)*sin(pi*r) - Vin*sin((r*pi*(2*Dp + Dy1 - Dy2 - 2))/2)*cos(pi*r))/(2*Z0*(cos(pi*r) + 1));
        % I_zvs_p34_RegionG = ((sin(pi*r)*(Vin + Vo*cos((r*pi*(2*Dp - Dy1 + Dy2 - 2))/2) - Vin*cos(pi*Dy1*r) - Vo*cos((r*pi*(Dy1 - 2*Dp + Dy2 + 2))/2)))/Z0 + ((cos(pi*r) + 1)*(Vo*sin((r*pi*(2*Dp - Dy1 + Dy2 - 2))/2) + Vin*sin(pi*Dy1*r) + Vo*sin((r*pi*(Dy1 - 2*Dp + Dy2 + 2))/2)))/Z0)/(2*(cos(pi*r) + 1));
        % I_zvs_s34_RegionG = -(((cos(pi*r) + 1)*(Vin*sin((r*pi*(2*Dp - Dy1 + Dy2))/2) - Vin*sin(pi*r) + Vo*sin(pi*Dy2*r) + Vin*sin((r*pi*(2*Dp + Dy1 + Dy2 - 2))/2)))/Z0 + (sin(pi*r)*(Vin + Vo - Vin*cos((r*pi*(2*Dp - Dy1 + Dy2))/2) + Vin*cos(pi*r) - Vo*cos(pi*Dy2*r) - Vin*cos((r*pi*(2*Dp + Dy1 + Dy2 - 2))/2)))/Z0)/(2*(cos(pi*r) + 1));
        % 
        % I_zvs_p12_RegionH = -(Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2) + Vin.*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2) + Vin.*sin(pi.*r.*(Dy1 - 1)) + Vo.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2).*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2).*cos(pi.*r) + Vin.*cos(pi.*r.*(Dy1 - 1)).*sin(pi.*r) + Vin.*sin(pi.*r.*(Dy1 - 1)).*cos(pi.*r) + Vo.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2).*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2).*cos(pi.*r))./(2.*Z0.*(cos(pi.*r) + 1));
        % I_zvs_s12_RegionH = (Vo.*sin(pi.*r) + Vo.*sin(pi.*Dy2.*r) + Vo.*cos(pi.*r).*sin(pi.*Dy2.*r) - Vo.*sin(pi.*r).*cos(pi.*Dy2.*r) + 2.*Vin.*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vin.*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vin.*cos(pi.*r).*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vin.*cos(pi.*r).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) - 2.*Vin.*sin(pi.*r).*cos(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vin.*sin(pi.*r).*cos((pi.*Dy2.*r)/2).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2))./(2.*Z0.*(cos(pi.*r) + 1));
        % I_zvs_p34_RegionH = ((sin(pi.*r).*(Vin + Vo.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) - Vin.*cos(pi.*Dy1.*r) - Vo.*cos((r.*pi.*(Dy1 - 2.*Dp + Dy2 + 2))/2)))/Z0 + ((cos(pi.*r) + 1).*(Vo.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) + Vin.*sin(pi.*Dy1.*r) + Vo.*sin((r.*pi.*(Dy1 - 2.*Dp + Dy2 + 2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));
        % I_zvs_s34_RegionH = -(((cos(pi.*r) + 1).*(Vin.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2))/2) - Vin.*sin(pi.*r) + Vo.*sin(pi.*Dy2.*r) + Vin.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2)))/Z0 + (sin(pi.*r).*(Vin + Vo - Vin.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2))/2) + Vin.*cos(pi.*r) - Vo.*cos(pi.*Dy2.*r) - Vin.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));

        Imax1 = abs((Vin.^2.*cos(pi.*r.*(Dy1 - 1)) + Vo.^2.*cos(pi.*r.*(Dy2 - 1)) + Vin.^2 + Vo.^2 - Vin.*Vo.*cos((pi.*r.*(2.*Dp + Dy1 + Dy2 - 2))/2) - Vin.*Vo.*cos((pi.*r.*(2.*Dp + Dy1 - Dy2))/2) - Vin.*Vo.*cos((pi.*r.*(2.*Dp - Dy1 + Dy2))/2) - Vin.*Vo.*cos((pi.*r.*(2.*Dp - Dy1 - Dy2 + 2))/2))./(Z0.^2.*(cos(pi.*r) + 1)));
        Imax2 = abs((2.*Vin.^2 - 2.*Vo.^2.*cos((Dy2.*pi.*r)/2).^2 - 2.*Vin.^2.*cos((Dy1.*pi.*r)/2).^2 + 2.*Vo.^2 + 4.*Vin.*Vo.*cos(pi.*r).*cos(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2) + 4.*Vin.*Vo.*sin(pi.*r).*sin(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2))./(Z0.^2.*(cos(pi.*r) + 1)));

    % calulate objective function
    switch region_flag
        case 1
            % Region A
            Ip12 = abs(I_zvs_p12_RegionA)/Vin;
            Is12 = abs(I_zvs_s12_RegionA)/Vin;  
            Ip34 = abs(I_zvs_p34_RegionA)/Vin;
            Is34 = abs(I_zvs_s34_RegionA)/Vin;
            ILmax = sqrt(Imax1)/Vin;
        case 2
            % Region B
            Ip12 = abs(I_zvs_p12_RegionB)/Vin;
            Is12 = abs(I_zvs_s12_RegionB)/Vin;
            Ip34 = abs(I_zvs_p34_RegionB)/Vin;
            Is34 = abs(I_zvs_s34_RegionB)/Vin;
            ILmax = sqrt(Imax1)/Vin;
        case 3
            % Region C
            Ip12 = abs(I_zvs_p12_RegionC)/Vin;
            Is12 = abs(I_zvs_s12_RegionC)/Vin;
            Ip34 = abs(I_zvs_p34_RegionC)/Vin;
            Is34 = abs(I_zvs_s34_RegionC)/Vin;
            ILmax = sqrt(Imax1)/Vin;
        case 4
            % Region D
            Ip12 = abs(I_zvs_p12_RegionD)/Vin;
            Is12 = abs(I_zvs_s12_RegionD)/Vin;
            Ip34 = abs(I_zvs_p34_RegionD)/Vin;
            Is34 = abs(I_zvs_s34_RegionD)/Vin;
            ILmax = sqrt(Imax1)/Vin;
        case 5
            % Region E
            Ip12 = abs(I_zvs_p12_RegionE)/Vin;
            Is12 = abs(I_zvs_s12_RegionE)/Vin;
            Ip34 = abs(I_zvs_p34_RegionE)/Vin;
            Is34 = abs(I_zvs_s34_RegionE)/Vin;
            ILmax = sqrt(Imax1)/Vin;
        % case 6
        %     % Region F
        %     h = -(((cos(pi.*r) + 1).*(Vo.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2 - 4))/2) + Vo.*sin(pi.*r) - Vin.*sin(pi.*Dy1.*r) - Vo.*sin((r.*pi.*(Dy1 - 2.*Dp + Dy2 + 2))/2)))/Z0 - (sin(pi.*r).*(Vin + Vo - Vo.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2 - 4))/2) + Vo.*cos(pi.*r) - Vin.*cos(pi.*Dy1.*r) - Vo.*cos((r.*pi.*(Dy1 - 2.*Dp + Dy2 + 2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));
        % case 7
        %     % Region G
        %     h = ((sin(pi*r)*(Vin + Vo*cos((r*pi*(2*Dp - Dy1 + Dy2 - 2))/2) - Vin*cos(pi*Dy1*r) - Vo*cos((r*pi*(Dy1 - 2*Dp + Dy2 + 2))/2)))/Z0 + ((cos(pi*r) + 1)*(Vo*sin((r*pi*(2*Dp - Dy1 + Dy2 - 2))/2) + Vin*sin(pi*Dy1*r) + Vo*sin((r*pi*(Dy1 - 2*Dp + Dy2 + 2))/2)))/Z0)/(2*(cos(pi*r) + 1));
        % case 8
        %     % Region H
        %     h = ((sin(pi.*r).*(Vin + Vo.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) - Vin.*cos(pi.*Dy1.*r) - Vo.*cos((r.*pi.*(Dy1 - 2.*Dp + Dy2 + 2))/2)))/Z0 + ((cos(pi.*r) + 1).*(Vo.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) + Vin.*sin(pi.*Dy1.*r) + Vo.*sin((r.*pi.*(Dy1 - 2.*Dp + Dy2 + 2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));
        otherwise
            error('Unknown region');
    end    
end


% 在parfor循环结束后，检测并替换突变点
function [optimal_points, optimal_power, optimal_Irms, Ioffp12, Ioffp34, Ioffs12, Ioffs34, ILmax, actual_PoN_values] = ...
    remove_outlier_points(optimal_points, optimal_power, optimal_Irms, Ioffp12, Ioffp34, Ioffs12, Ioffs34, ILmax, actual_PoN_values, Vin, Vo_values, DoF)
    
    % 首先处理无效点（值为0的点）
    valid_indices = ~any(optimal_points == 0, 2);
    
    % 然后检测突变点
    [m, ~, p] = size(optimal_points);
    outlier_mask = false(size(optimal_points, 1), size(optimal_points, 3));
    
    % 定义突变阈值（根据实际数据调整）
    threshold_dp = 0.05;
    threshold_dy = 0.1;
    threshold_r = 0.005;
    
    % 对每个Vo值检测
    for v_idx = 1:p
        % 对每个PoN值检测突变
        for i = 2:m-1  % 跳过第一个和最后一个点
            % 检查当前点与相邻点的差异
            diff_prev = abs(optimal_points(i,:,v_idx) - optimal_points(i-1,:,v_idx));
            diff_next = abs(optimal_points(i,:,v_idx) - optimal_points(i+1,:,v_idx));
            
            % 应用不同阈值到不同参数
            thresholds = [threshold_dp, threshold_dy, threshold_dy];
            if DoF >= 4
                thresholds = [thresholds, threshold_r];
            end
            
            % 如果差异过大，标记为突变点
            if any(diff_prev > thresholds) && any(diff_next > thresholds)
                outlier_mask(i, v_idx) = true;
            end
        end
        
        % 处理突变点 - 使用相邻点的平均值替换
        for i = 2:m-1
            if outlier_mask(i, v_idx)
                % 查找前后有效的非突变点
                prev_valid = max(1, i-1);
                while prev_valid > 1 && (outlier_mask(prev_valid, v_idx) || ~valid_indices(prev_valid, v_idx))
                    prev_valid = prev_valid - 1;
                end
                
                next_valid = min(m, i+1);
                while next_valid < m && (outlier_mask(next_valid, v_idx) || ~valid_indices(next_valid, v_idx))
                    next_valid = next_valid + 1;
                end
                
                % 确保找到了有效点
                if valid_indices(prev_valid, v_idx) && valid_indices(next_valid, v_idx)
                    % 线性插值替换突变点
                    weight = (i - prev_valid) / (next_valid - prev_valid);
                    interp_vars = (1-weight) * optimal_points(prev_valid,:,v_idx) + weight * optimal_points(next_valid,:,v_idx);
                    
                    % 更新所有相关值
                    optimal_points(i,:,v_idx) = interp_vars;
                    Vo = Vo_values(v_idx);
                    
                    % 使用插值参数重新计算物理量
                    optimal_power(i,v_idx) = f_function(DoF, interp_vars, Vin, Vo);
                    optimal_Irms(i,v_idx) = g_function(DoF, interp_vars, Vin, Vo);
                    [Ip12_val, Ip34_val, Is12_val, Is34_val, ILmax_val] = h_function(DoF, interp_vars, Vin, Vo);
                    Ioffp12(i,v_idx) = Ip34_val;
                    Ioffp34(i,v_idx) = Ip12_val;
                    Ioffs12(i,v_idx) = Is34_val;
                    Ioffs34(i,v_idx) = Is12_val;
                    ILmax(i,v_idx) = ILmax_val;
                    actual_PoN_values(i,v_idx) = f_function(DoF, interp_vars, Vin, Vo);
                end
            end
        end
    end
end