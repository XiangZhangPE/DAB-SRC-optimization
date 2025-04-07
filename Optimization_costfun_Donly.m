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
Vo = 660;  % Changed from 1000 to match g_function definition

plotoptions = 2;
ZVSPlot_DABsrc(r,Vo,Vin,Z0, plotoptions); hold on;

% Define optimization problem
options = optimoptions('fmincon', 'Display', 'off');
% options = optimoptions('fmincon', 'Display', 'final');
% options = optimoptions('fmincon', ...
%     'OutputFcn', @g_function, ...  % Output during iteration
%     'Display', 'iter', ...          % Print each iteration in the command window
%     'ConstraintTolerance', 1e-9, ...
%     'OptimalityTolerance', 1e-6, ... 
%     'StepTolerance', 1e-9, ...
%     'MaxIterations', 100);          % Adjustable maximum iterations



% bounds = [0 0 0; 1 1 0.5];      % Lower and upper bounds
% initial_guess = [0.5, 0.5, 0.25];  % Initial guess
% PoN_values = linspace(0.5, 1.75, 100);

bounds = [0 0 0; 1 1 0.5];      % Lower and upper bounds
initial_guess = [0.5, 0.5, 0.25];  % Initial guess
PoN_values = linspace(0.45, 1.75, 100);

%% Run optimization for each PoN target value
optimal_points = zeros(length(PoN_values), 3);
optimal_power = zeros(length(PoN_values), 1);
optimal_Irms = zeros(length(PoN_values), 1);

for i = 1:length(PoN_values)
    PoN_target = PoN_values(i);
    
    % Define objective function with penalty
    obj_fun = @(vars) J_function(vars, r, Vin, Vo, Z0, PoN_target);  % Objective function now includes PoN_target
    
    % Define nonlinear constraints, import r, Vin, Vo, Z0
    nonlcon = @(vars) zvs_constrain(vars, r, Vin, Vo, Z0);
    
    % Run optimization
    [optimal_vars, fval, exitflag, output] = fmincon(obj_fun, initial_guess, [], [], [], [], bounds(1,:), bounds(2,:), nonlcon, options);

    % Store result
    optimal_points(i,:) = optimal_vars;
    optimal_power(i,:) = f_function(optimal_vars, r);
    optimal_Irms(i,:) = g_function(optimal_vars, r, Vin, Vo, Z0);
    
    % Update initial guess for next optimization
    initial_guess = optimal_vars;
    
    % Output optimization status
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
view([2 5 3]); % M<1
% view([5 2 3]); % M>1

figure
plot(PoN_values, optimal_Irms, 'b-', 'LineWidth', 1); hold on;
title('Irms under different Power','Interpreter', 'latex');
xlabel('Normalized output power PoN', 'Interpreter', 'latex');
ylabel('Normalized tank current Irms', 'Interpreter', 'latex');
grid on;

figure
plot(PoN_values, optimal_points(:,1), 'k-', 'LineWidth', 1, 'DisplayName', 'Dy1'); hold on;
plot(PoN_values, optimal_points(:,2), 'r-', 'LineWidth', 1, 'DisplayName', 'Dy2'); hold on;
plot(PoN_values, optimal_points(:,3), 'b-', 'LineWidth', 1, 'DisplayName', 'Dyp'); hold on;
xlabel('Normalized output power PoN', 'Interpreter', 'latex');
ylabel('Phase-shift ratios','Interpreter', 'latex');
legend('show');
grid on;


% Define Cost function with penalty for PoN_target
function J = J_function(vars, r, Vin, Vo, Z0, PoN_target)

    PoN = f_function(vars, r);
    Irms = g_function(vars, r, Vin, Vo, Z0);

    % Add penalty term (adjustable weight λ)
    lambda = 500;  % Controls penalty strength
    penalty = lambda * (PoN - PoN_target)^2;
    
    % Final cost function
    J = Irms + penalty;
end

% Define g(Dy1, Dy2, Dp) function with penalty for PoN_target
function g = g_function(vars, r, Vin, Vo, Z0)
    Dy1 = vars(1);
    Dy2 = vars(2);
    Dp = vars(3);
    
    % Original objective function calculation
    g = Z0/Vin*((1/2/pi/r)*(-(2*Vin^2*sin(Dy1*pi*r) + 2*Vo^2*sin(Dy2*pi*r) + 2*Vin^2*sin(pi*r*(Dy1 - 1)) + 2*Vo^2*sin(pi*r*(Dy2 - 1)) + 2*Vin^2*sin(pi*r) + 2*Vo^2*sin(pi*r) + 2*Vin*Vo*sin((pi*r*(2*Dp - Dy1 - Dy2 + 2))/2) - 2*Vin*Vo*sin((pi*r*(Dy1 - 2*Dp + Dy2))/2) + 2*Vin*Vo*sin((pi*r*(2*Dp + Dy1 + Dy2 - 2))/2) + 2*Vin*Vo*sin((pi*r*(2*Dp + Dy1 + Dy2 - 4))/2) + 2*Vin*Vo*sin((pi*r*(2*Dp + Dy1 - Dy2))/2) + 2*Vin*Vo*sin((pi*r*(2*Dp - Dy1 + Dy2))/2) + 2*Vin*Vo*sin((pi*r*(2*Dp + Dy1 - Dy2 - 2))/2) + 2*Vin*Vo*sin((pi*r*(2*Dp - Dy1 + Dy2 - 2))/2) - 2*Vin^2*pi*r - 2*Vo^2*pi*r + 2*Vin^2*pi*r*cos(Dy1*pi*r) + 2*Vo^2*pi*r*cos(Dy2*pi*r) - 2*Vin*Vo*pi*r*cos((pi*r*(Dy1 - 2*Dp + Dy2))/2) + 4*Vin*Vo*pi*r*cos((pi*r*(2*Dp + Dy1 + Dy2 - 2))/2) + 2*Vin*Vo*pi*r*cos((pi*r*(2*Dp + Dy1 + Dy2 - 4))/2) - 2*Dy1*Vin^2*pi*r*cos(Dy1*pi*r) - 2*Dy2*Vo^2*pi*r*cos(Dy2*pi*r) + 2*Vin*Vo*pi*r*cos((pi*r*(2*Dp + Dy1 - Dy2))/2) + 2*Vin*Vo*pi*r*cos((pi*r*(2*Dp - Dy1 + Dy2))/2) - 2*Dy1*Vin^2*pi*r*cos(pi*r*(Dy1 - 1)) - 2*Dy2*Vo^2*pi*r*cos(pi*r*(Dy2 - 1)) - 2*Dp*Vin*Vo*pi*r*cos((pi*r*(2*Dp - Dy1 - Dy2 + 2))/2) + Dy1*Vin*Vo*pi*r*cos((pi*r*(2*Dp - Dy1 - Dy2 + 2))/2) + Dy2*Vin*Vo*pi*r*cos((pi*r*(2*Dp - Dy1 - Dy2 + 2))/2) - 2*Dp*Vin*Vo*pi*r*cos((pi*r*(Dy1 - 2*Dp + Dy2))/2) - 2*Dp*Vin*Vo*pi*r*cos((pi*r*(2*Dp + Dy1 + Dy2 - 2))/2) - 2*Dp*Vin*Vo*pi*r*cos((pi*r*(2*Dp + Dy1 + Dy2 - 4))/2) + Dy1*Vin*Vo*pi*r*cos((pi*r*(Dy1 - 2*Dp + Dy2))/2) + Dy2*Vin*Vo*pi*r*cos((pi*r*(Dy1 - 2*Dp + Dy2))/2) - Dy1*Vin*Vo*pi*r*cos((pi*r*(2*Dp + Dy1 + Dy2 - 2))/2) - Dy2*Vin*Vo*pi*r*cos((pi*r*(2*Dp + Dy1 + Dy2 - 2))/2) - Dy1*Vin*Vo*pi*r*cos((pi*r*(2*Dp + Dy1 + Dy2 - 4))/2) - Dy2*Vin*Vo*pi*r*cos((pi*r*(2*Dp + Dy1 + Dy2 - 4))/2))/(4*Z0^2*(cos(pi*r) + 1))));
end

% Define f(Dy1, Dy2, Dp) function
function f = f_function(vars, r)
    Dy1 = vars(1);
    Dy2 = vars(2);
    Dp = vars(3);

    % Calculate f value
    f = ((cos(pi*r*(Dp-0.5))*cos(pi*r*(1-Dy1)/2)*cos(pi*r*(1-Dy2)/2))/(pi*r*cos(pi*r/2))-1/(pi*r));
end

% Define constraint function for ZVS
function [c, ceq] = zvs_constrain(vars, r, Vin, Vo, Z0)
    Dy1 = vars(1);
    Dy2 = vars(2);
    Dp = vars(3);
    
    % ZVS constraints
    I_zvs_p12_RegionD = -(Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2) + Vin.*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2) + Vin.*sin(pi.*r.*(Dy1 - 1)) + Vo.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2).*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2).*cos(pi.*r) + Vin.*cos(pi.*r.*(Dy1 - 1)).*sin(pi.*r) + Vin.*sin(pi.*r.*(Dy1 - 1)).*cos(pi.*r) + Vo.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2).*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2).*cos(pi.*r))./(2.*Z0.*(cos(pi.*r) + 1));
    I_zvs_s12_RegionD = (Vin.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2))/2) + Vo.*sin(pi.*r) + Vin.*sin((r.*pi.*(2.*Dp - Dy1 - Dy2 + 2))/2) + Vo.*sin(pi.*r.*(Dy2 - 1)) + Vin.*cos(pi.*r).*sin((r.*pi.*(2.*Dp - Dy1 - Dy2 + 2))/2) - Vin.*sin(pi.*r).*cos((r.*pi.*(2.*Dp - Dy1 - Dy2 + 2))/2) + Vo.*cos(pi.*r.*(Dy2 - 1)).*sin(pi.*r) + Vo.*sin(pi.*r.*(Dy2 - 1)).*cos(pi.*r) - Vin.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2))/2).*sin(pi.*r) + Vin.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2))/2).*cos(pi.*r))./(2.*Z0.*(cos(pi.*r) + 1));
    I_zvs_p34_RegionD =  (((cos(pi.*r) + 1).*(Vo.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) + Vo.*sin(pi.*r) + Vin.*sin(pi.*Dy1.*r) - Vo.*sin((r.*pi.*(Dy1 - 2.*Dp + Dy2))/2)))/Z0 + (sin(pi.*r).*(Vin - Vo + Vo.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) - Vo.*cos(pi.*r) - Vin.*cos(pi.*Dy1.*r) + Vo.*cos((r.*pi.*(Dy1 - 2.*Dp + Dy2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));
    I_zvs_s34_RegionD = -(((cos(pi.*r) + 1).*(Vin.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2))/2) - Vin.*sin(pi.*r) + Vo.*sin(pi.*Dy2.*r) + Vin.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2)))/Z0 + (sin(pi.*r).*(Vin + Vo - Vin.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2))/2) + Vin.*cos(pi.*r) - Vo.*cos(pi.*Dy2.*r) - Vin.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));
     
    % Define inequality constraints (c <= 0)
    c = [ I_zvs_p12_RegionD;    % I_zvs_p12_RegionD <= 0
         -I_zvs_s12_RegionD;    % -I_zvs_s12_RegionD <= 0 (equivalent to I_zvs_s12_RegionD >= 0)
         -I_zvs_p34_RegionD;    % -I_zvs_p34_RegionD <= 0 (equivalent to I_zvs_p34_RegionD >= 0)
         I_zvs_s34_RegionD];    % I_zvs_s34_RegionD <= 0
         
    % Define equality constraint (ceq = 0)
    ceq = [];     % PoN should equal PoN_target
end
