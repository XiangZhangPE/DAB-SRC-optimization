% Define constants
fs = 250e3;
Lr = 12e-6;
Cr = 48e-9;
Z0 = sqrt(Lr/Cr);
fr = 1/(2*pi*sqrt(Lr*Cr));
td = 50e-9;
Dd = td*fs;
r = fr/fs;
omega_r = 2*pi*fr;
Ln = 100;
Vin = 800;
Vo = 580;  % Changed from 1000 to match g_function definition

plotoptions = 2;
ZVSPlot_DABsrc(r,Vo,Vin,Z0, plotoptions); hold on;

% Define optimization problem
options = optimoptions('fmincon', 'Display', 'off');
bounds = [0 0 0; 1 1 0.5];      % Lower and upper bounds


% Run optimization for each PoN target value

% low high power guess
initial_guess = [0.7, 0.95, 0.06];  % Initial guess
PoN_values = linspace(0.015, 0.45, 100);

% % % meduim to high power guess
% initial_guess = [1, 1, 0.25];  % Initial guess
% PoN_values = linspace(0.5, 1.75, 100);

optimal_points = zeros(length(PoN_values), 3);
optimal_power = zeros(length(PoN_values), 3);
optimal_Irms = zeros(length(PoN_values), 3);

for i = 1:length(PoN_values)
    PoN_target = PoN_values(i);
    
    % Define objective function with penalty
    obj_fun = @(vars) J_function(vars, r, Vin, Vo, Z0, PoN_target);  % 目标函数现在包含 PoN_target
    
    % Define nonlinear constraints, 传入 r, Vin, Vo, Z0
    nonlcon = @(vars) zvs_constrain(vars, r, Vin, Vo, Z0);
    
    % Run optimization
    [optimal_vars, fval, exitflag, output] = fmincon(obj_fun, initial_guess, [], [], [], [], bounds(1,:), bounds(2,:), nonlcon, options);

    % Store result
    optimal_points(i,:) = optimal_vars;
    optimal_power(i,:) = f_function(optimal_vars, r);
    optimal_Irms(i,:) = g_function(optimal_vars, r, Vin, Vo, Z0);
    
    % Update initial guess for next optimization
    initial_guess = optimal_vars;
    
    % 输出优化状态
    if exitflag > 0
        fprintf('✅ Optimization successful for PoN_target = %.4f. Final Cost: %.6f\n', PoN_target, fval);
        fprintf('   Optimal Variables: Dy1 = %.6f, Dy2 = %.6f, Dp = %.6f\n', optimal_vars(1), optimal_vars(2), optimal_vars(3));
    else
        warning('⚠️ Optimization did NOT converge for PoN_target = %.4f! Exit flag: %d\n', PoN_target, exitflag);
        fprintf('   Reason: %s\n', output.message);
        fprintf('   Last Attempt Variables: Dy1 = %.6f, Dy2 = %.6f, Dp = %.6f\n', optimal_vars(1), optimal_vars(2), optimal_vars(3));
    end
end


% Visualize results in 3D
scatter3(optimal_points(:,1), optimal_points(:,2), optimal_points(:,3), 10, 'red','filled'); % 50 is the point size
hold on;
plot3(optimal_points(:,1), optimal_points(:,2), optimal_points(:,3), 'red', 'LineWidth', 2);
xlabel('Dy1');
ylabel('Dy2');
zlabel('Dp');
title('Optimized g(Dy1, Dy2, Dp) for different PoN values');
xlim([0 1]);
ylim([0 1]);
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

figure
plot(optimal_power, optimal_Irms);



% Define g(Dy1, Dy2, Dp) function with penalty for PoN_target
function J = J_function(vars, r, Vin, Vo, Z0, PoN_target)

    % 计算 PoN
    PoN = f_function(vars, r);
    % 计算Irms
    Irms = g_function(vars, r, Vin, Vo, Z0);
    
    % 添加惩罚项（权重 λ 可调整）
    lambda = 1000;  % 控制惩罚强度
    penalty = lambda * (PoN - PoN_target)^2;
    
    % 最终目标函数
    J = Irms + penalty;
end



