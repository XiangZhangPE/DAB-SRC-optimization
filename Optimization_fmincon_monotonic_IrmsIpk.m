% clc, close all;
%% Define constants
fs = 250e3; % switching frequency, user defined
Lr = 12e-6; % resonant inductor, user defined
Cr = 48e-9; % resonant capacitor, user defined
Z0 = sqrt(Lr/Cr);
fr = 1/(2*pi*sqrt(Lr*Cr));
r = fr/fs;

Vin = 100; % input voltage, user defined
Vo = 80;  % output voltage, user defined

% unused variables
omega_r = 2*pi*fr;
td = 50e-9;
Dd = td*fs;
Ln = 100;
ILm0 = pi*r*Vo/(2 * Ln * Z0);
Ioss = 0.06;

% plot background ZVS figure
handle = figure(3);
plotoptions = 2;
ZVSPlot_DABsrc(r,Vo,Vin,Z0, plotoptions); hold on;

%% Optimization problem
% Define the power ranges 
% low_power_range = linspace(0.1, 0.45, 50); % power range for Ipk optimization
low_power_range = linspace(0.015, 0.45, 50);
high_power_range = linspace(0.45, 1.2, 50);
PoN_values = [low_power_range, high_power_range]; % Include both power ranges
DoF = 4;

% Define different initial conditions for each region
high_power_initial = [0.15, 1, 1, 1/1.1];
low_power_buck_initial = [0.012, 0.29, 0.38, 0.77];
low_power_boost_initial = [0.012, 0.38, 0.29, 0.77];
if (Vo>=Vin)
    low_power_initial = low_power_boost_initial;
else
    low_power_initial = low_power_buck_initial;
end

% Define bounds for optimization variables based on DoF
if DoF == 1
    bounds = [0; 0.5];  % Lower and upper bounds for Dp
elseif DoF == 2
    bounds = [0 0; 0.5 1];  % Bounds for Dp, Dy1
elseif DoF == 3
    bounds = [0 0 0; 0.5 1 1];  % Bounds for Dp, Dy1, Dy2
elseif DoF == 4
    % bounds = [0 0 0 1/1.4; 0.5 1 1 1/1.01];  % Bounds for Ipk optimization
    bounds = [0 0 0 1/1.3; 0.5 1 1 1/1.05];  % Bounds for Dp, Dy1, Dy2, r
end

% Prepare arrays with extra space
total_points = length(PoN_values);
optimal_points = zeros(total_points, DoF);
optimal_power = zeros(total_points, 1);
optimal_Irms = zeros(total_points, 1);
Ioffp12 = zeros(total_points, 1);
Ioffp34 = zeros(total_points, 1);
Ioffs12 = zeros(total_points, 1);
Ioffs34 = zeros(total_points, 1);
ILmax = zeros(total_points, 1);
actual_PoN_values = zeros(total_points, 1);

% Update optimization options to increase iterations
options = optimoptions('fmincon', ...
    'Display', 'off', ...
    'MaxIterations', 2000, ...  % Increase maximum iterations
    'MaxFunctionEvaluations', 10000, ...  % Increase function evaluations
    'OptimalityTolerance', 1e-8, ...
    'StepTolerance', 1e-10);

% define threshold for two adjacent optimal variables
distance_threshold = 0.02; % maximum allowed steps between two optimals
distance_threshold1 = -0.001; % maximum decremental step
distance_threshold2 = 0.02; % maximum incremental step

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
            fval = fval_low;
            exitflag = exitflag_low;
            fprintf('✅ Optimization successful for PoN_target = %.4f (low power region selected)\n', PoN_target);
        else
            optimal_vars = optimal_vars_high;
            fval = fval_high;
            exitflag = exitflag_high;
            fprintf('✅ Optimization successful for PoN_target = %.4f (high power region selected)\n', PoN_target);
        end
    elseif exitflag_low > 0
        optimal_vars = optimal_vars_low;
        fval = fval_low;
        exitflag = exitflag_low;
        fprintf('✅ Optimization successful for PoN_target = %.4f (low power region only)\n', PoN_target);
    elseif exitflag_high > 0
        optimal_vars = optimal_vars_high;
        fval = fval_high;
        exitflag = exitflag_high;
        fprintf('✅ Optimization successful for PoN_target = %.4f (high power region only)\n', PoN_target);
    else
        warning('⚠️ Both optimizations did NOT converge for PoN_target = %.4f', PoN_target);
        continue; % Skip to next iteration
    end

    % % 检查点是否减少过多 (允许减少但有下限)
    % if i > 1
    %     % 计算当前点与前一个点的差值
    %     diff_vars = optimal_vars - optimal_points(i-1, :);
        
    %     % 限制减少的下限
    %     excessive_decrement = diff_vars < distance_threshold1;
    %     optimal_vars(excessive_decrement) = optimal_points(i-1, excessive_decrement) + distance_threshold1;
        
    %     % 限制增量不超过阈值
    %     excessive_increment = diff_vars > distance_threshold2;
    %     optimal_vars(excessive_increment) = optimal_points(i-1, excessive_increment) + distance_threshold2;

    %     % 打印调整信息
    %     if any(excessive_decrement)
    %         fprintf('⚠️ Excessive increment detected at PoN_target = %.4f. Limited.\n', PoN_target);
    %     end
    % end
    
    % checkpoint for monotonic optimization
    % adjacent not over incremental step limit
    % matrix based method
    if i > 1
        % calculate the adjacent variables' difference
        diff_vars = optimal_vars - optimal_points(i-1, :);
        
        % force monotonic increasing
        optimal_vars(diff_vars < 0) = optimal_points(i-1, diff_vars < 0);
        
        % incremental step limits for adjacent variables
        excessive_increment = diff_vars > distance_threshold;
        optimal_vars(excessive_increment) = optimal_points(i-1, excessive_increment) + distance_threshold;
        
        % print warning if violate the rules
        if any(diff_vars < 0)
            fprintf('⚠️ Non-monotonic behavior detected at PoN_target = %.4f. Adjusted to maintain monotonicity.\n', PoN_target);
        end
        if any(excessive_increment)
            fprintf('⚠️ Excessive increment detected at PoN_target = %.4f. Increment limited.\n', PoN_target);
        end
    end
    % update inital guess for next iteration
    low_power_initial = optimal_vars; % use last successful optimization result for next iteration
    high_power_initial = optimal_vars; % use last successful optimization result for next iteration

    
    % Store results
    optimal_points(i,:) = optimal_vars;
    optimal_power(i,:) = f_function(DoF, optimal_vars, Vin, Vo);
    optimal_Irms(i,:) = g_function(DoF, optimal_vars, Vin, Vo);
    [Ip12_val, Ip34_val, Is12_val, Is34_val, ILmax_val] = h_function(DoF, optimal_vars, Vin, Vo);
    Ioffp12(i,:) = Ip34_val;
    Ioffp34(i,:) = Ip12_val;
    Ioffs12(i,:) = Is34_val;
    Ioffs34(i,:) = Is12_val;
    ILmax(i,:) = ILmax_val;
    actual_PoN_values(i) = f_function(DoF, optimal_vars, Vin, Vo); % 使用 f_function 计算实际 PoN 值
end

% Remove any uninitialized points (where optimization failed)
valid_indices = actual_PoN_values ~= 0;
optimal_points = optimal_points(valid_indices, :);
optimal_power = optimal_power(valid_indices);
optimal_Irms = optimal_Irms(valid_indices);
Ioffp12 = Ioffp12(valid_indices);
Ioffp34 = Ioffp34(valid_indices);
Ioffs12 = Ioffs12(valid_indices);
Ioffs34 = Ioffs34(valid_indices);
ILmax = ILmax(valid_indices);
actual_PoN_values = actual_PoN_values(valid_indices);

% Use actual PoN values for plotting
figure(1)
plot(actual_PoN_values, Ioffp12, 'r-', 'LineWidth', 1, 'DisplayName', 'I_{off,P12}'); hold on;
plot(actual_PoN_values, Ioffp34, 'b-', 'LineWidth', 1, 'DisplayName', 'I_{off,P34}'); hold on;
plot(actual_PoN_values, Ioffs12, 'm-', 'LineWidth', 1, 'DisplayName', 'I_{off,S12}'); hold on;
plot(actual_PoN_values, Ioffs34, 'c-', 'LineWidth', 1, 'DisplayName', 'I_{off,S34}'); hold on;
plot(actual_PoN_values, ILmax, 'k--', 'LineWidth', 1, 'DisplayName', 'I_{L,max}'); hold on;
ylabel('Normalized turn-off current', 'Interpreter', 'latex');
title('Ioff under different Power','Interpreter', 'latex');
xlabel('Normalized output power PoN', 'Interpreter', 'latex');
legend('show');
grid on;

figure(2)
plot(actual_PoN_values, optimal_points(:,1), 'k-', 'LineWidth', 1, 'DisplayName', 'Outer-phase-shift Dp'); hold on;
plot(actual_PoN_values, optimal_points(:,2), 'r-', 'LineWidth', 1, 'DisplayName', 'Inner-phase-shift Dy1'); hold on;
plot(actual_PoN_values, optimal_points(:,3), 'b-', 'LineWidth', 1, 'DisplayName', 'Inner-phase-shift Dy2'); hold on;
if DoF >= 4
    plot(actual_PoN_values, optimal_points(:,4), 'g-', 'LineWidth', 1, 'DisplayName', 'Frequency ratio r = fr/fs'); hold on;
end
xlabel('Normalized output power PoN', 'Interpreter', 'latex');
ylabel('Phase-shift ratios','Interpreter', 'latex');
legend('show');
grid on;

figure(3)
if DoF >= 3
    % Visualize results in 3D
    scatter3(optimal_points(:,2), optimal_points(:,3), optimal_points(:,1), 10, 'red','filled');
    hold on;
    plot3(optimal_points(:,2), optimal_points(:,3), optimal_points(:,1), 'red', 'LineWidth', 2);
    xlabel('Dy1');
    ylabel('Dy2');
    zlabel('Dp');
    title('Optimized g(Dy1, Dy2, Dp) for different PoN values');
    xlim([0 1]);
    ylim([0 1]);
    % zlim([0 1]);
    zlim([0 0.5]);
    grid on;
    if (Vo/Vin<1)
        view([2 5 3]); % M<1
    elseif (Vo/Vin>1)
        view([5 2 3]); % M>1
    else
        view([4 4 5]);
    end
    % view([-160, 40]);
elseif DoF == 2
    % Visualize results in 2D
    scatter(optimal_points(:,2), optimal_points(:,1), 10, 'red','filled'); % 50 is the point size
    hold on;
    plot(optimal_points(:,2), optimal_points(:,1), 'red', 'LineWidth', 2);
    xlabel('Dy1');
    ylabel('Dp');
    title('Optimized g(Dy1, Dp) for different PoN values');
    xlim([0 1]);
    ylim([0 0.5]);
    grid on;
end





% J_function
function J = J_function(DoF, vars, Vin, Vo, PoN_target)
    % calculate PoN
    PoN = f_function(DoF, vars, Vin, Vo);
    % calculate Irms
    Irms = g_function(DoF, vars, Vin, Vo);
    Ipk = g2_function(DoF, vars, Vin, Vo);
    % calculate Ioff (you can ignore the return value if not needed)
    h_function(DoF, vars, Vin, Vo);
    
    % construct cost function
    lambda = 800;  % add penalty coefficient (weight λ adjustable)
    penalty = lambda * (PoN - PoN_target)^2;
    
    % final cost function
    % J = Irms + penalty;
    J = Ipk + penalty;
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
    region_flag((((Dy2+Dy1)/2)<= Dp) & (Dp<(1-(Dy2+Dy1)/2))) = 5;  %  E
    region_flag(((1-(Dy2+Dy1)/2)<= Dp) & (((Dy2+Dy1)/2)>= Dp)) = 4;  %  D
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
        case 5
            % Region E
            f = 1/2.*(((sin(pi*r.*Dy1/2).*cos(pi*r.*(2*Dp-1)/2).*sin(pi*r.*Dy2/2))./(r.*cos(pi*r/2))));
        case 4
            % Region D
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


% Define Ipk function
function g2 = g2_function(DoF, vars, Vin, Vo)
    Z0 = 1;
    % Get variables based on DoF
    [Dp, Dy1, Dy2, r] = assign_variables(DoF, vars, Vin, Vo);

    region_flag((Dy1>Dy2) & (Dp<(Dy1-Dy2)/2)) = 1;  %  A
    region_flag((Dy2>=Dy1) & (Dp<(Dy2-Dy1)/2)) = 2;  %  B
    region_flag((Dp>=abs(Dy1-Dy2)/2) & (Dp<(Dy2+Dy1)/2) & (Dp<(1-(Dy2+Dy1)/2))) = 3;  %  C
    region_flag((((Dy2+Dy1)/2)<= Dp) & (Dp<(1-(Dy2+Dy1)/2))) = 5;  %  E
    region_flag(((1-(Dy2+Dy1)/2)<= Dp) & (((Dy2+Dy1)/2)>= Dp)) = 4;  %  D
    % region_flag(((1-(Dy2-Dy1)/2)<=Dp) & (Dy2>Dy1)) = 6;  %  F
    % region_flag(((1-(Dy1-Dy2)/2)<=Dp) & (Dy2<=Dy1)) = 7;  %  G
    % region_flag((Dp>=(1-(Dy1+Dy2)/2)) & (Dp>=(Dy1+Dy2)/2) & (Dp<(1-abs(Dy1-Dy2)/2))) = 8;  %  H

    % calulate objective function
    if region_flag<=4
        Imax_2 = (2*abs(Vo.^2.*(cos(pi.*r.*(1-Dy2)/2)).^2+Vin.^2.*(cos(pi.*r.*(1-Dy1)/2)).^2-2.*Vin.*Vo.*cos(pi.*r.*(Dp)).*cos(pi.*r.*(1-Dy1)/2).*cos(pi.*r.*(1-Dy2)/2)))./(Z0.^2.*(cos(pi.*r) + 1));
    else
        Imax_2 = (2*abs(Vo.^2.*(sin(pi.*r.*(Dy2)/2)).^2+Vin.^2.*(sin(pi.*r.*(Dy1)/2)).^2 + 2.*Vin.*Vo.*cos(pi.*r.*(1+Dp)).*sin(pi.*r.*(Dy1)/2).*sin(pi.*r.*(Dy2)/2)))./(Z0.^2.*(cos(pi.*r) + 1));
    end
    g2 = sqrt(Imax_2)/Vin;
end


% Define Irms function
function g = g_function(DoF, vars, Vin, Vo)
    % Get variables based on DoF
    [Dp, Dy1, Dy2, r] = assign_variables(DoF, vars, Vin, Vo);

    region_flag((Dy1>Dy2) & (Dp<(Dy1-Dy2)/2)) = 1;  %  A
    region_flag((Dy2>=Dy1) & (Dp<(Dy2-Dy1)/2)) = 2;  %  B
    region_flag((Dp>=abs(Dy1-Dy2)/2) & (Dp<(Dy2+Dy1)/2) & (Dp<(1-(Dy2+Dy1)/2))) = 3;  %  C
    region_flag((((Dy2+Dy1)/2)<= Dp) & (Dp<(1-(Dy2+Dy1)/2))) = 5;  %  E
    region_flag(((1-(Dy2+Dy1)/2)<= Dp) & (((Dy2+Dy1)/2)>= Dp)) = 4;  %  D
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
        case 5
            % Region E
            g = 0;
        case 4
            % Region D
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
    region_flag((((Dy2+Dy1)/2)<= Dp) & (Dp<(1-(Dy2+Dy1)/2))) = 5;  %  E
    region_flag(((1-(Dy2+Dy1)/2)<= Dp) & (((Dy2+Dy1)/2)>= Dp)) = 4;  %  D
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

    I_zvs_p12_RegionE = -((sin(pi.*r).*(Vo.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2) - Vin.*cos(pi.*r) - Vo.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2) + Vin.*cos(pi.*r.*(Dy1 - 1))))/Z0 + ((cos(pi.*r) + 1).*(Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2) + Vin.*sin(pi.*r) - Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2) + Vin.*sin(pi.*r.*(Dy1 - 1))))/Z0)./(2.*(cos(pi.*r) + 1));
    I_zvs_s12_RegionE = (Vo.*sin(pi.*r) + Vo.*sin(pi.*Dy2.*r) + Vo.*cos(pi.*r).*sin(pi.*Dy2.*r) - Vo.*sin(pi.*r).*cos(pi.*Dy2.*r) + 2.*Vin.*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vin.*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vin.*cos(pi.*r).*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vin.*cos(pi.*r).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) - 2.*Vin.*sin(pi.*r).*cos(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vin.*sin(pi.*r).*cos((pi.*Dy2.*r)/2).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2))./(2.*Z0.*(cos(pi.*r) + 1));
    I_zvs_p34_RegionE = ((sin(pi.*r).*(Vin + Vo.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) - Vin.*cos(pi.*Dy1.*r) - Vo.*cos((r.*pi.*(Dy1 - 2.*Dp + Dy2 + 2))/2)))/Z0 + ((cos(pi.*r) + 1).*(Vo.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) + Vin.*sin(pi.*Dy1.*r) + Vo.*sin((r.*pi.*(Dy1 - 2.*Dp + Dy2 + 2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));
    I_zvs_s34_RegionE = ((sin(pi.*r).*(2.*Vo.*cos((pi.*Dy2.*r)/2).^2 - 2.*Vo + 2.*Vin.*cos(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vin.*cos((pi.*Dy2.*r)/2).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2)))/Z0 - ((cos(pi.*r) + 1).*(2.*Vo.*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy2.*r)/2) - 2.*Vin.*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vin.*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2)))/Z0)./(2.*(cos(pi.*r) + 1));

    I_zvs_p12_RegionD = -(Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2) + Vin.*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2) + Vin.*sin(pi.*r.*(Dy1 - 1)) + Vo.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2).*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2).*cos(pi.*r) + Vin.*cos(pi.*r.*(Dy1 - 1)).*sin(pi.*r) + Vin.*sin(pi.*r.*(Dy1 - 1)).*cos(pi.*r) + Vo.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2).*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2).*cos(pi.*r))./(2.*Z0.*(cos(pi.*r) + 1));
    I_zvs_s12_RegionD = (Vin.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2))/2) + Vo.*sin(pi.*r) + Vin.*sin((r.*pi.*(2.*Dp - Dy1 - Dy2 + 2))/2) + Vo.*sin(pi.*r.*(Dy2 - 1)) + Vin.*cos(pi.*r).*sin((r.*pi.*(2.*Dp - Dy1 - Dy2 + 2))/2) - Vin.*sin(pi.*r).*cos((r.*pi.*(2.*Dp - Dy1 - Dy2 + 2))/2) + Vo.*cos(pi.*r.*(Dy2 - 1)).*sin(pi.*r) + Vo.*sin(pi.*r.*(Dy2 - 1)).*cos(pi.*r) - Vin.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2))/2).*sin(pi.*r) + Vin.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2))/2).*cos(pi.*r))./(2.*Z0.*(cos(pi.*r) + 1));
    I_zvs_p34_RegionD =  (((cos(pi.*r) + 1).*(Vo.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) + Vo.*sin(pi.*r) + Vin.*sin(pi.*Dy1.*r) - Vo.*sin((r.*pi.*(Dy1 - 2.*Dp + Dy2))/2)))/Z0 + (sin(pi.*r).*(Vin - Vo + Vo.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) - Vo.*cos(pi.*r) - Vin.*cos(pi.*Dy1.*r) + Vo.*cos((r.*pi.*(Dy1 - 2.*Dp + Dy2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));
    I_zvs_s34_RegionD = -(((cos(pi.*r) + 1).*(Vin.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2))/2) - Vin.*sin(pi.*r) + Vo.*sin(pi.*Dy2.*r) + Vin.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2)))/Z0 + (sin(pi.*r).*(Vin + Vo - Vin.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2))/2) + Vin.*cos(pi.*r) - Vo.*cos(pi.*Dy2.*r) - Vin.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));

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
    % 
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
        case 5
            % Region E
            c = [I_zvs_p12_RegionD; -I_zvs_s12_RegionD; -I_zvs_p34_RegionD; I_zvs_s34_RegionD];
        case 4
            % Region D
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
    region_flag((((Dy2+Dy1)/2)<= Dp) & (Dp<(1-(Dy2+Dy1)/2))) = 5;  %  E
    region_flag(((1-(Dy2+Dy1)/2)<= Dp) & (((Dy2+Dy1)/2)>= Dp)) = 4;  %  D
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
    
        I_zvs_p12_RegionE = -((sin(pi.*r).*(Vo.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2) - Vin.*cos(pi.*r) - Vo.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2) + Vin.*cos(pi.*r.*(Dy1 - 1))))/Z0 + ((cos(pi.*r) + 1).*(Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2) + Vin.*sin(pi.*r) - Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2) + Vin.*sin(pi.*r.*(Dy1 - 1))))/Z0)./(2.*(cos(pi.*r) + 1));
        I_zvs_s12_RegionE = (Vo.*sin(pi.*r) + Vo.*sin(pi.*Dy2.*r) + Vo.*cos(pi.*r).*sin(pi.*Dy2.*r) - Vo.*sin(pi.*r).*cos(pi.*Dy2.*r) + 2.*Vin.*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vin.*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vin.*cos(pi.*r).*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vin.*cos(pi.*r).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) - 2.*Vin.*sin(pi.*r).*cos(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vin.*sin(pi.*r).*cos((pi.*Dy2.*r)/2).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2))./(2.*Z0.*(cos(pi.*r) + 1));
        I_zvs_p34_RegionE = ((sin(pi.*r).*(Vin + Vo.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) - Vin.*cos(pi.*Dy1.*r) - Vo.*cos((r.*pi.*(Dy1 - 2.*Dp + Dy2 + 2))/2)))/Z0 + ((cos(pi.*r) + 1).*(Vo.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) + Vin.*sin(pi.*Dy1.*r) + Vo.*sin((r.*pi.*(Dy1 - 2.*Dp + Dy2 + 2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));
        I_zvs_s34_RegionE = ((sin(pi.*r).*(2.*Vo.*cos((pi.*Dy2.*r)/2).^2 - 2.*Vo + 2.*Vin.*cos(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vin.*cos((pi.*Dy2.*r)/2).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2)))/Z0 - ((cos(pi.*r) + 1).*(2.*Vo.*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy2.*r)/2) - 2.*Vin.*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vin.*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2)))/Z0)./(2.*(cos(pi.*r) + 1));
    
        I_zvs_p12_RegionD = -(Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2) + Vin.*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2) + Vin.*sin(pi.*r.*(Dy1 - 1)) + Vo.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2).*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2).*cos(pi.*r) + Vin.*cos(pi.*r.*(Dy1 - 1)).*sin(pi.*r) + Vin.*sin(pi.*r.*(Dy1 - 1)).*cos(pi.*r) + Vo.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2).*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2).*cos(pi.*r))./(2.*Z0.*(cos(pi.*r) + 1));
        I_zvs_s12_RegionD = (Vin.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2))/2) + Vo.*sin(pi.*r) + Vin.*sin((r.*pi.*(2.*Dp - Dy1 - Dy2 + 2))/2) + Vo.*sin(pi.*r.*(Dy2 - 1)) + Vin.*cos(pi.*r).*sin((r.*pi.*(2.*Dp - Dy1 - Dy2 + 2))/2) - Vin.*sin(pi.*r).*cos((r.*pi.*(2.*Dp - Dy1 - Dy2 + 2))/2) + Vo.*cos(pi.*r.*(Dy2 - 1)).*sin(pi.*r) + Vo.*sin(pi.*r.*(Dy2 - 1)).*cos(pi.*r) - Vin.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2))/2).*sin(pi.*r) + Vin.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2))/2).*cos(pi.*r))./(2.*Z0.*(cos(pi.*r) + 1));
        I_zvs_p34_RegionD =  (((cos(pi.*r) + 1).*(Vo.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) + Vo.*sin(pi.*r) + Vin.*sin(pi.*Dy1.*r) - Vo.*sin((r.*pi.*(Dy1 - 2.*Dp + Dy2))/2)))/Z0 + (sin(pi.*r).*(Vin - Vo + Vo.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) - Vo.*cos(pi.*r) - Vin.*cos(pi.*Dy1.*r) + Vo.*cos((r.*pi.*(Dy1 - 2.*Dp + Dy2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));
        I_zvs_s34_RegionD = -(((cos(pi.*r) + 1).*(Vin.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2))/2) - Vin.*sin(pi.*r) + Vo.*sin(pi.*Dy2.*r) + Vin.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2)))/Z0 + (sin(pi.*r).*(Vin + Vo - Vin.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2))/2) + Vin.*cos(pi.*r) - Vo.*cos(pi.*Dy2.*r) - Vin.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));
    
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

        % Imax1 = abs((Vin.^2.*cos(pi.*r.*(Dy1 - 1)) + Vo.^2.*cos(pi.*r.*(Dy2 - 1)) + Vin.^2 + Vo.^2 - Vin.*Vo.*cos((pi.*r.*(2.*Dp + Dy1 + Dy2 - 2))/2) - Vin.*Vo.*cos((pi.*r.*(2.*Dp + Dy1 - Dy2))/2) - Vin.*Vo.*cos((pi.*r.*(2.*Dp - Dy1 + Dy2))/2) - Vin.*Vo.*cos((pi.*r.*(2.*Dp - Dy1 - Dy2 + 2))/2))./(Z0.^2.*(cos(pi.*r) + 1)));
        % Imax2 = abs((2.*Vin.^2 - 2.*Vo.^2.*cos((Dy2.*pi.*r)/2).^2 - 2.*Vin.^2.*cos((Dy1.*pi.*r)/2).^2 + 2.*Vo.^2 + 4.*Vin.*Vo.*cos(pi.*r).*cos(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2) + 4.*Vin.*Vo.*sin(pi.*r).*sin(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2))./(Z0.^2.*(cos(pi.*r) + 1)));
        
        Imax1 = (2*abs(Vo.^2.*(cos(pi.*r.*(1-Dy2)/2)).^2+Vin.^2.*(cos(pi.*r.*(1-Dy1)/2)).^2-2.*Vin.*Vo.*cos(pi.*r.*(Dp)).*cos(pi.*r.*(1-Dy1)/2).*cos(pi.*r.*(1-Dy2)/2)))./(Z0.^2.*(cos(pi.*r) + 1));
        Imax2 = (2*abs(Vo.^2.*(sin(pi.*r.*(Dy2)/2)).^2+Vin.^2.*(sin(pi.*r.*(Dy1)/2)).^2 + 2.*Vin.*Vo.*cos(pi.*r.*(1+Dp)).*sin(pi.*r.*(Dy1)/2).*sin(pi.*r.*(Dy2)/2)))./(Z0.^2.*(cos(pi.*r) + 1));

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
            ILmax = sqrt(Imax2)/Vin;
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