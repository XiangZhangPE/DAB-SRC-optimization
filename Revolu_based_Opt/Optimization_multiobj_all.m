clc, close all;
%% Define constants - unchanged from original code
fs = 250e3; % switching frequency, user defined
Lr = 12e-6; % resonant inductor, user defined
Cr = 48e-9; % resonant capacitor, user defined
Z0 = sqrt(Lr/Cr);
fr = 1/(2*pi*sqrt(Lr*Cr));
r = fr/fs;

Vin = 800; % input voltage, user defined
Vo = 1200;  % output voltage, user defined

% unused variables
omega_r = 2*pi*fr;
td = 50e-9;
Dd = td*fs;
Ln = 100;
ILm0 = pi*r*Vo/(2 * Ln * Z0);
Ioss = 0.06;

% Plot background ZVS figure (optional)
handle = figure(3);
plotoptions = 2;
ZVSPlot_DABsrc(r,Vo,Vin,Z0, plotoptions); hold on;

%% Setup for multi-objective optimization
% Define bounds and degrees of freedom

% bounds = [0 0 0; 0.5 1 1]; % Lower and upper bounds for DoF=3
% DoF = 3 % Degree of freedom: set to 3

bounds = [0 0 0 1/1.2; 0.5 1 1 1/1.1];
DoF = 4;

% Create initial population for genetic algorithm with appropriate dimensionality
popSize = 200;

% Create base population uniformly distributed within bounds
initPop = repmat(bounds(1,:), popSize/2, 1) + ...
          rand(popSize/2, DoF) .* repmat((bounds(2,:) - bounds(1,:)), popSize/2, 1);

% Add targeted initial points based on DoF
switch DoF
    case 1
        % Only Dp is variable
        targetPoints1 = [0.25] + 0.1*randn(popSize/4, 1);
        targetPoints2 = [0.4] + 0.05*randn(popSize/4, 1);
        
    case 2
        % Dp and one of Dy1/Dy2 are variable
        if Vin > Vo
            targetPoints1 = [0.25, 0.8] + 0.1*randn(popSize/4, 2);
            targetPoints2 = [0.4, 0.9] + 0.05*randn(popSize/4, 2);
        else
            targetPoints1 = [0.25, 0.8] + 0.1*randn(popSize/4, 2);
            targetPoints2 = [0.4, 0.9] + 0.05*randn(popSize/4, 2);
        end
        
    case 3
        % Dp, Dy1, and Dy2 are variable
        targetPoints1 = [0.25, 0.8, 0.8] + 0.1*randn(popSize/4, 3);
        targetPoints2 = [0.4, 1, 1] + 0.05*randn(popSize/4, 3);
        
    case 4
        % All four parameters are variable
        targetPoints1 = [0.25, 0.8, 0.8, 1/1.15] + 0.1*randn(popSize/4, 4);
        targetPoints2 = [0.4, 1, 1, 1/1.2] + 0.05*randn(popSize/4, 4);
end

% Combine all initial points
initPop = [initPop; targetPoints1; targetPoints2];

% Set options for multi-objective optimization
% options = optimoptions('gamultiobj', ...
%     'PopulationSize', popSize, ...
%     'InitialPopulation', initPop, ...
%     'Display', 'iter', ...
%     'MaxGenerations', 200, ...
%     'ParetoFraction', 0.7, ...
%     'UseParallel', true);  % Enable parallel processing if available
options = optimoptions('gamultiobj', ...
    'PopulationSize', 200,...            % 增加种群规模
    'InitialPopulation', initPop, ...
    'Display', 'iter', ...
    'MaxGenerations', 800, ...                    % 增加最大迭代次数
    'ParetoFraction', 0.4, ...                    % 调整Pareto保留比例
    'CrossoverFraction', 0.8, ...                 % 调整交叉概率
    'MutationFcn', @mutationadaptfeasible, ...    % 使用自适应变异函数
    'CrossoverFcn', @crossoverintermediate,  ...  % 尝试不同的交叉函数
    'FunctionTolerance', 1e-5, ...                % 更严格的函数收敛容差
    'UseParallel', true);

%% Run multi-objective optimization
% Define multi-objective function directly invoking existing g and f functions
multiObj_fun = @(vars) [g_function(DoF, vars, Vin, Vo), 1./f_function(DoF, vars, Vin, Vo)];

% Define nonlinear constraints, unchanged
nonlcon = @(vars) zvs_function(DoF, vars, Vin, Vo, Z0);

% Run optimization to find Pareto front
[optimal_vars, objective_values] = gamultiobj(multiObj_fun, DoF, [], [], [], [], ...
    bounds(1,:), bounds(2,:), nonlcon, options);

%% Process and visualize results
% Extract Irms and PoN values from objectives
Irms_values = objective_values(:,1);
PoN_values = 1./objective_values(:,2);

% Calculate turn-off currents for the Pareto-optimal solutions
num_solutions = size(optimal_vars, 1);
Ioffp12 = zeros(num_solutions, 1);
Ioffp34 = zeros(num_solutions, 1);
Ioffs12 = zeros(num_solutions, 1);
Ioffs34 = zeros(num_solutions, 1);
ILmax = zeros(num_solutions, 1);

for i = 1:num_solutions
    [Ip12, Ip34, Is12, Is34, IL] = h_function(DoF, optimal_vars(i,:), Vin, Vo);
    Ioffp12(i) = Ip12;
    Ioffp34(i) = Ip34;
    Ioffs12(i) = Is12;
    Ioffs34(i) = Is34;
    ILmax(i) = IL;
end

%% Visualization of Pareto front
figure(1);
scatter(PoN_values, Irms_values, 5, 'filled');
xlabel('Normalized output power PoN', 'Interpreter', 'latex');
ylabel('Normalized tank current I_{rms}', 'Interpreter', 'latex');
title('Pareto Front: I_{rms} vs. PoN', 'Interpreter', 'latex');
grid on;

% Sort solutions by PoN for line plotting
[PoN_sorted, idx] = sort(PoN_values);
Irms_sorted = Irms_values(idx);
optimal_vars_sorted = optimal_vars(idx,:);

figure(2);
plot(PoN_sorted, optimal_vars_sorted(:,1), 'k-', 'LineWidth', 1, 'DisplayName', 'Dp'); hold on;
plot(PoN_sorted, optimal_vars_sorted(:,2), 'r-', 'LineWidth', 1, 'DisplayName', 'Dy1'); hold on;
plot(PoN_sorted, optimal_vars_sorted(:,3), 'b-', 'LineWidth', 1, 'DisplayName', 'Dy2'); hold on;
if DoF >= 4
    plot(PoN_sorted, optimal_vars_sorted(:,4), 'g-', 'LineWidth', 1, 'DisplayName', 'r'); hold on;
end
xlabel('Normalized output power PoN', 'Interpreter', 'latex');
ylabel('Phase-shift ratios', 'Interpreter', 'latex');
legend('show');
grid on;

figure(3);
% Visualize the Pareto optimal points on the ZVS plot
if DoF >= 3
    scatter3(optimal_vars(:,2), optimal_vars(:,3), optimal_vars(:,1), 5, 'r', 'filled');
    colorbar;
    xlabel('Dy1');
    ylabel('Dy2');
    zlabel('Dp');
    title('Pareto Optimal Solutions on ZVS Plot');   
    view([45, 30]);
elseif DoF == 2
    scatter(optimal_vars(:,2), optimal_vars(:,1), 5, 'r', 'filled');
    colorbar;
    xlabel('Dy1');
    ylabel('Dp');
    title('Pareto Optimal Solutions on ZVS Plot');
end

% Turn-off current visualization
figure(4);
plot(PoN_sorted, Ioffp12(idx), 'r-', 'LineWidth', 1, 'DisplayName', 'I_{off,P12}'); hold on;
plot(PoN_sorted, Ioffp34(idx), 'b-', 'LineWidth', 1, 'DisplayName', 'I_{off,P34}'); hold on;
plot(PoN_sorted, Ioffs12(idx), 'm-', 'LineWidth', 1, 'DisplayName', 'I_{off,S12}'); hold on;
plot(PoN_sorted, Ioffs34(idx), 'c-', 'LineWidth', 1, 'DisplayName', 'I_{off,S34}'); hold on;
plot(PoN_sorted, ILmax(idx), 'k--', 'LineWidth', 1, 'DisplayName', 'I_{L,max}'); hold on;
xlabel('Normalized output power PoN', 'Interpreter', 'latex');
ylabel('Normalized turn-off current', 'Interpreter', 'latex');
title('Turn-off Currents Along Pareto Front', 'Interpreter', 'latex');
legend('show');
grid on;

%% Analysis of Pareto front
% Find knee point (using "elbow" method)
% Normalize objectives to 0-1 range
Irms_norm = (Irms_values - min(Irms_values)) / (max(Irms_values) - min(Irms_values));
PoN_norm = (PoN_values - min(PoN_values)) / (max(PoN_values) - min(PoN_values));

% Calculate distances from origin to each point
distances = sqrt(Irms_norm.^2 + (1-PoN_norm).^2);
[~, knee_idx] = min(distances);

% Highlight knee point
figure(1); hold on;
scatter(PoN_values(knee_idx), Irms_values(knee_idx), 100, 'r', 'filled', 'DisplayName', 'Knee Point');
legend('Pareto Front', 'Knee Point');

% Report results based on DoF
fprintf('Knee point solution:\n');
fprintf('PoN = %.4f, Irms = %.4f\n', PoN_values(knee_idx), Irms_values(knee_idx));

switch DoF
    case 1
        fprintf('Dp = %.4f\n', optimal_vars(knee_idx,1));
    case 2
        fprintf('Dp = %.4f, Dy1/Dy2 = %.4f\n', optimal_vars(knee_idx,1), optimal_vars(knee_idx,2));
    case 3
        fprintf('Dp = %.4f, Dy1 = %.4f, Dy2 = %.4f\n', optimal_vars(knee_idx,1), optimal_vars(knee_idx,2), optimal_vars(knee_idx,3));
    case 4
        fprintf('Dp = %.4f, Dy1 = %.4f, Dy2 = %.4f, r = %.4f\n', optimal_vars(knee_idx,1), optimal_vars(knee_idx,2), optimal_vars(knee_idx,3), optimal_vars(knee_idx,4));
end

% Save results, including DoF for reference
save('pareto_results.mat', 'optimal_vars', 'Irms_values', 'PoN_values', 'knee_idx', 'DoF');






% Define power function
function f = f_function(DoF, vars, Vin, Vo)

    % Dp = vars(1); Dy1 = vars(2); Dy2 = vars(3); r = 1/1.2;
    if DoF >= 1
        Dp = vars(1); Dy1 = 1; Dy2 = 1; r = 1/1.2;
    end
    if DoF >= 2
        if (Vin>Vo)
            Dy1 = vars(2);
            Dy2 = 1;
        else
            Dy1 = 1;
            Dy2 = vars(2);
        end
        r = 1/1.2;
    end
    if DoF >= 3
        Dy1 = vars(2); Dy2 = vars(3); r = 1/1.2;
    end
    if DoF >= 4
        r = vars(4);
    end

    region_flag((Dy1>Dy2) & (Dp<(Dy1-Dy2)/2)) = 1;  %  A
    region_flag((Dy2>=Dy1) & (Dp<(Dy2-Dy1)/2)) = 2;  %  B
    region_flag((Dp>=abs(Dy1-Dy2)/2) & (Dp<(Dy2+Dy1)/2) & (Dp<(1-(Dy2+Dy1)/2))) = 3;  %  C
    region_flag((((Dy2+Dy1)/2)<= Dp) & (Dp<(1-(Dy2+Dy1)/2))) = 4;  %  D
    region_flag(((1-(Dy2+Dy1)/2)<= Dp) & (((Dy2+Dy1)/2)>= Dp)) = 5;  %  E
    % region_flag(((1-(Dy2-Dy1)/2)<=Dp) & (Dy2>Dy1)) = 6;  %  F
    % region_flag(((1-(Dy1-Dy2)/2)<=Dp) & (Dy2<=Dy1)) = 7;  %  G
    % region_flag((Dp>=(1-(Dy1+Dy2)/2)) & (Dp>=(Dy1+Dy2)/2) & (Dp<(1-abs(Dy1-Dy2)/2))) = 8;  %  H

    % Calculate f value
    % switch region_flag
    %     case 1
    %         % Region A
    %         f =  1/pi.*((sin(pi*r.*Dy2/2).*sin(pi*r.*Dp).*cos(pi*r.*(1-Dy1)/2))./(r.*cos(pi*r/2)));
    %     case 2
    %         % Region B
    %         f = 1/pi.*((sin(pi*r.*Dy1/2).*sin(pi*r.*Dp).*cos(pi*r.*(1-Dy2)/2))./(r.*cos(pi*r/2)));
    %     case 3
    %         % Region C
    %         f = 1/(2*pi).*(((cos(pi.*r.*(Dy2-Dy1)/2).*cos(pi.*r/2.*(2.*Dp-1)) + sin(pi.*r.*(Dy1+Dy2-1)/2).*sin(pi.*r.*(Dp)))./(r.*cos(pi.*r./2)) - 1./r));
    %     case 4
    %         % Region D
    %         f = 1/pi.*(((sin(pi*r.*Dy1/2).*cos(pi*r.*(2*Dp-1)/2).*sin(pi*r.*Dy2/2))./(r.*cos(pi*r/2))));
    %     case 5
    %         % Region E
    %         f = 1/pi.*((cos(pi.*r.*(Dp-0.5)).*cos(pi.*r.*(1-Dy1)./2).*cos(pi.*r.*(1-Dy2)./2))./(r.*cos(pi.*r./2))-1./(r));
    %     % case 6
    %     %     % Region F
    %     %     f = 1/pi.*((sin(pi*r.*Dy1/2).*sin(pi*r.*(1-Dp)).*cos(pi*r.*(1-Dy2)/2))./(r.*cos(pi*r/2)));
    %     % case 7
    %     %     % Region G
    %     %     f = 1/pi.*((sin(pi*r.*Dy2/2).*sin(pi*r.*(1-Dp)).*cos(pi*r.*(1-Dy1)/2))./(r.*cos(pi*r/2)));
    %     % case 8
    %     %     % Region H
    %     %     f = 1/(2*pi).*(((cos(pi.*r.*(Dy2-Dy1)/2).*cos(pi.*r/2.*(2.*Dp-1)) + sin(pi.*r.*(Dy1+Dy2-1)/2).*sin(pi.*r.*(1-Dp)))./(r.*cos(pi.*r./2)) - 1./r));
    %     otherwise
    %         error('Unknown region');
    % end

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

    % Dp = vars(1); Dy1 = vars(2); Dy2 = vars(3); r = 1/1.2;
    if DoF >= 1
        Dp = vars(1); Dy1 = 1; Dy2 = 1; r = 1/1.2;
    end
    if DoF >= 2
        if (Vin>Vo)
            Dy1 = vars(2);
            Dy2 = 1;
        else
            Dy1 = 1;
            Dy2 = vars(2);
        end
        r = 1/1.2;
    end
    if DoF >= 3
        Dy1 = vars(2); Dy2 = vars(3); r = 1/1.2;
    end
    if DoF >= 4
        r = vars(4);
    end

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

    % Dp = vars(1); Dy1 = vars(2); Dy2 = vars(3); r = 1/1.2;
    if DoF >= 1
        Dp = vars(1); Dy1 = 1; Dy2 = 1; r = 1/1.2;
    end
    if DoF >= 2
        if (Vin>Vo)
            Dy1 = vars(2);
            Dy2 = 1;
        else
            Dy1 = 1;
            Dy2 = vars(2);
        end
        r = 1/1.2;
    end
    if DoF >= 3
        Dy1 = vars(2); Dy2 = vars(3); r = 1/1.2;
    end
    if DoF >= 4
        r = vars(4);
    end

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

    I_zvs_p12_RegionF = -(Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2) + Vin.*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2) + Vin.*sin(pi.*r.*(Dy1 - 1)) + Vo.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2).*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2).*cos(pi.*r) + Vin.*cos(pi.*r.*(Dy1 - 1)).*sin(pi.*r) + Vin.*sin(pi.*r.*(Dy1 - 1)).*cos(pi.*r) + Vo.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2).*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2).*cos(pi.*r))./(2.*Z0.*(cos(pi.*r) + 1));
    I_zvs_s12_RegionF = (Vo.*sin(pi.*r) + Vo.*sin(pi.*Dy2.*r) + Vo.*cos(pi.*r).*sin(pi.*Dy2.*r) - Vo.*sin(pi.*r).*cos(pi.*Dy2.*r) + 2.*Vin.*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vin.*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vin.*cos(pi.*r).*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vin.*cos(pi.*r).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) - 2.*Vin.*sin(pi.*r).*cos(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vin.*sin(pi.*r).*cos((pi.*Dy2.*r)/2).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2))./(2.*Z0.*(cos(pi.*r) + 1));
    I_zvs_p34_RegionF = -(((cos(pi.*r) + 1).*(Vo.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2 - 4))/2) + Vo.*sin(pi.*r) - Vin.*sin(pi.*Dy1.*r) - Vo.*sin((r.*pi.*(Dy1 - 2.*Dp + Dy2 + 2))/2)))/Z0 - (sin(pi.*r).*(Vin + Vo - Vo.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2 - 4))/2) + Vo.*cos(pi.*r) - Vin.*cos(pi.*Dy1.*r) - Vo.*cos((r.*pi.*(Dy1 - 2.*Dp + Dy2 + 2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));
    I_zvs_s34_RegionF = -((sin(pi.*r).*(Vo + Vin.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) - Vo.*cos(pi.*Dy2.*r) - Vin.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2)))/Z0 + ((cos(pi.*r) + 1).*(Vo.*sin(pi.*Dy2.*r) - Vin.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) + Vin.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));

    I_zvs_p12_RegionG = ((sin(pi*r)*(Vo*cos((r*pi*(2*Dp + Dy1 - Dy2 - 4))/2) + Vin*cos(pi*r) - Vo*cos((r*pi*(2*Dp + Dy1 + Dy2 - 4))/2) - Vin*cos(pi*r*(Dy1 - 1))))/Z0 - ((cos(pi*r) + 1)*(Vin*sin(pi*r) - Vo*sin((r*pi*(2*Dp + Dy1 - Dy2 - 4))/2) + Vo*sin((r*pi*(2*Dp + Dy1 + Dy2 - 4))/2) + Vin*sin(pi*r*(Dy1 - 1))))/Z0)/(2*(cos(pi*r) + 1));
    I_zvs_s12_RegionG = (Vo*sin(pi*r) - Vin*sin((r*pi*(2*Dp + Dy1 - Dy2 - 2))/2) + Vin*sin((r*pi*(Dy1 - 2*Dp + Dy2))/2) + Vo*sin(pi*r*(Dy2 - 1)) + Vin*cos((r*pi*(Dy1 - 2*Dp + Dy2))/2)*sin(pi*r) + Vin*sin((r*pi*(Dy1 - 2*Dp + Dy2))/2)*cos(pi*r) + Vo*cos(pi*r*(Dy2 - 1))*sin(pi*r) + Vo*sin(pi*r*(Dy2 - 1))*cos(pi*r) + Vin*cos((r*pi*(2*Dp + Dy1 - Dy2 - 2))/2)*sin(pi*r) - Vin*sin((r*pi*(2*Dp + Dy1 - Dy2 - 2))/2)*cos(pi*r))/(2*Z0*(cos(pi*r) + 1));
    I_zvs_p34_RegionG = ((sin(pi*r)*(Vin + Vo*cos((r*pi*(2*Dp - Dy1 + Dy2 - 2))/2) - Vin*cos(pi*Dy1*r) - Vo*cos((r*pi*(Dy1 - 2*Dp + Dy2 + 2))/2)))/Z0 + ((cos(pi*r) + 1)*(Vo*sin((r*pi*(2*Dp - Dy1 + Dy2 - 2))/2) + Vin*sin(pi*Dy1*r) + Vo*sin((r*pi*(Dy1 - 2*Dp + Dy2 + 2))/2)))/Z0)/(2*(cos(pi*r) + 1));
    I_zvs_s34_RegionG = -(((cos(pi*r) + 1)*(Vin*sin((r*pi*(2*Dp - Dy1 + Dy2))/2) - Vin*sin(pi*r) + Vo*sin(pi*Dy2*r) + Vin*sin((r*pi*(2*Dp + Dy1 + Dy2 - 2))/2)))/Z0 + (sin(pi*r)*(Vin + Vo - Vin*cos((r*pi*(2*Dp - Dy1 + Dy2))/2) + Vin*cos(pi*r) - Vo*cos(pi*Dy2*r) - Vin*cos((r*pi*(2*Dp + Dy1 + Dy2 - 2))/2)))/Z0)/(2*(cos(pi*r) + 1));

    I_zvs_p12_RegionH = -(Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2) + Vin.*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2) + Vin.*sin(pi.*r.*(Dy1 - 1)) + Vo.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2).*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2).*cos(pi.*r) + Vin.*cos(pi.*r.*(Dy1 - 1)).*sin(pi.*r) + Vin.*sin(pi.*r.*(Dy1 - 1)).*cos(pi.*r) + Vo.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2).*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2).*cos(pi.*r))./(2.*Z0.*(cos(pi.*r) + 1));
    I_zvs_s12_RegionH = (Vo.*sin(pi.*r) + Vo.*sin(pi.*Dy2.*r) + Vo.*cos(pi.*r).*sin(pi.*Dy2.*r) - Vo.*sin(pi.*r).*cos(pi.*Dy2.*r) + 2.*Vin.*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vin.*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vin.*cos(pi.*r).*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vin.*cos(pi.*r).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) - 2.*Vin.*sin(pi.*r).*cos(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vin.*sin(pi.*r).*cos((pi.*Dy2.*r)/2).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2))./(2.*Z0.*(cos(pi.*r) + 1));
    I_zvs_p34_RegionH = ((sin(pi.*r).*(Vin + Vo.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) - Vin.*cos(pi.*Dy1.*r) - Vo.*cos((r.*pi.*(Dy1 - 2.*Dp + Dy2 + 2))/2)))/Z0 + ((cos(pi.*r) + 1).*(Vo.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) + Vin.*sin(pi.*Dy1.*r) + Vo.*sin((r.*pi.*(Dy1 - 2.*Dp + Dy2 + 2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));
    I_zvs_s34_RegionH = -(((cos(pi.*r) + 1).*(Vin.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2))/2) - Vin.*sin(pi.*r) + Vo.*sin(pi.*Dy2.*r) + Vin.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2)))/Z0 + (sin(pi.*r).*(Vin + Vo - Vin.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2))/2) + Vin.*cos(pi.*r) - Vo.*cos(pi.*Dy2.*r) - Vin.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));

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
    % Dp = vars(1); Dy1 = vars(2); Dy2 = vars(3); r = 1/1.2;
    if DoF >= 1
        Dp = vars(1); Dy1 = 1; Dy2 = 1; r = 1/1.2;
    end
    if DoF >= 2
        if (Vin>Vo)
            Dy1 = vars(2);
            Dy2 = 1;
        else
            Dy1 = 1;
            Dy2 = vars(2);
        end
        r = 1/1.2;
    end
    if DoF >= 3
        Dy1 = vars(2); Dy2 = vars(3); r = 1/1.2;
    end
    if DoF >= 4
        r = vars(4);
    end

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
    
        I_zvs_p12_RegionF = -(Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2) + Vin.*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2) + Vin.*sin(pi.*r.*(Dy1 - 1)) + Vo.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2).*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2).*cos(pi.*r) + Vin.*cos(pi.*r.*(Dy1 - 1)).*sin(pi.*r) + Vin.*sin(pi.*r.*(Dy1 - 1)).*cos(pi.*r) + Vo.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2).*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2).*cos(pi.*r))./(2.*Z0.*(cos(pi.*r) + 1));
        I_zvs_s12_RegionF = (Vo.*sin(pi.*r) + Vo.*sin(pi.*Dy2.*r) + Vo.*cos(pi.*r).*sin(pi.*Dy2.*r) - Vo.*sin(pi.*r).*cos(pi.*Dy2.*r) + 2.*Vin.*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vin.*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vin.*cos(pi.*r).*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vin.*cos(pi.*r).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) - 2.*Vin.*sin(pi.*r).*cos(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vin.*sin(pi.*r).*cos((pi.*Dy2.*r)/2).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2))./(2.*Z0.*(cos(pi.*r) + 1));
        I_zvs_p34_RegionF = -(((cos(pi.*r) + 1).*(Vo.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2 - 4))/2) + Vo.*sin(pi.*r) - Vin.*sin(pi.*Dy1.*r) - Vo.*sin((r.*pi.*(Dy1 - 2.*Dp + Dy2 + 2))/2)))/Z0 - (sin(pi.*r).*(Vin + Vo - Vo.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2 - 4))/2) + Vo.*cos(pi.*r) - Vin.*cos(pi.*Dy1.*r) - Vo.*cos((r.*pi.*(Dy1 - 2.*Dp + Dy2 + 2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));
        I_zvs_s34_RegionF = -((sin(pi.*r).*(Vo + Vin.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) - Vo.*cos(pi.*Dy2.*r) - Vin.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2)))/Z0 + ((cos(pi.*r) + 1).*(Vo.*sin(pi.*Dy2.*r) - Vin.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) + Vin.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));
    
        I_zvs_p12_RegionG = ((sin(pi*r)*(Vo*cos((r*pi*(2*Dp + Dy1 - Dy2 - 4))/2) + Vin*cos(pi*r) - Vo*cos((r*pi*(2*Dp + Dy1 + Dy2 - 4))/2) - Vin*cos(pi*r*(Dy1 - 1))))/Z0 - ((cos(pi*r) + 1)*(Vin*sin(pi*r) - Vo*sin((r*pi*(2*Dp + Dy1 - Dy2 - 4))/2) + Vo*sin((r*pi*(2*Dp + Dy1 + Dy2 - 4))/2) + Vin*sin(pi*r*(Dy1 - 1))))/Z0)/(2*(cos(pi*r) + 1));
        I_zvs_s12_RegionG = (Vo*sin(pi*r) - Vin*sin((r*pi*(2*Dp + Dy1 - Dy2 - 2))/2) + Vin*sin((r*pi*(Dy1 - 2*Dp + Dy2))/2) + Vo*sin(pi*r*(Dy2 - 1)) + Vin*cos((r*pi*(Dy1 - 2*Dp + Dy2))/2)*sin(pi*r) + Vin*sin((r*pi*(Dy1 - 2*Dp + Dy2))/2)*cos(pi*r) + Vo*cos(pi*r*(Dy2 - 1))*sin(pi*r) + Vo*sin(pi*r*(Dy2 - 1))*cos(pi*r) + Vin*cos((r*pi*(2*Dp + Dy1 - Dy2 - 2))/2)*sin(pi*r) - Vin*sin((r*pi*(2*Dp + Dy1 - Dy2 - 2))/2)*cos(pi*r))/(2*Z0*(cos(pi*r) + 1));
        I_zvs_p34_RegionG = ((sin(pi*r)*(Vin + Vo*cos((r*pi*(2*Dp - Dy1 + Dy2 - 2))/2) - Vin*cos(pi*Dy1*r) - Vo*cos((r*pi*(Dy1 - 2*Dp + Dy2 + 2))/2)))/Z0 + ((cos(pi*r) + 1)*(Vo*sin((r*pi*(2*Dp - Dy1 + Dy2 - 2))/2) + Vin*sin(pi*Dy1*r) + Vo*sin((r*pi*(Dy1 - 2*Dp + Dy2 + 2))/2)))/Z0)/(2*(cos(pi*r) + 1));
        I_zvs_s34_RegionG = -(((cos(pi*r) + 1)*(Vin*sin((r*pi*(2*Dp - Dy1 + Dy2))/2) - Vin*sin(pi*r) + Vo*sin(pi*Dy2*r) + Vin*sin((r*pi*(2*Dp + Dy1 + Dy2 - 2))/2)))/Z0 + (sin(pi*r)*(Vin + Vo - Vin*cos((r*pi*(2*Dp - Dy1 + Dy2))/2) + Vin*cos(pi*r) - Vo*cos(pi*Dy2*r) - Vin*cos((r*pi*(2*Dp + Dy1 + Dy2 - 2))/2)))/Z0)/(2*(cos(pi*r) + 1));
    
        I_zvs_p12_RegionH = -(Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2) + Vin.*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2) + Vin.*sin(pi.*r.*(Dy1 - 1)) + Vo.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2).*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2).*cos(pi.*r) + Vin.*cos(pi.*r.*(Dy1 - 1)).*sin(pi.*r) + Vin.*sin(pi.*r.*(Dy1 - 1)).*cos(pi.*r) + Vo.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2).*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2).*cos(pi.*r))./(2.*Z0.*(cos(pi.*r) + 1));
        I_zvs_s12_RegionH = (Vo.*sin(pi.*r) + Vo.*sin(pi.*Dy2.*r) + Vo.*cos(pi.*r).*sin(pi.*Dy2.*r) - Vo.*sin(pi.*r).*cos(pi.*Dy2.*r) + 2.*Vin.*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vin.*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vin.*cos(pi.*r).*cos(pi.*Dp.*r).*cos((pi.*Dy2.*r)/2).*sin((pi.*Dy1.*r)/2) + 2.*Vin.*cos(pi.*r).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) - 2.*Vin.*sin(pi.*r).*cos(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2).*sin((pi.*Dy2.*r)/2) + 2.*Vin.*sin(pi.*r).*cos((pi.*Dy2.*r)/2).*sin(pi.*Dp.*r).*sin((pi.*Dy1.*r)/2))./(2.*Z0.*(cos(pi.*r) + 1));
        I_zvs_p34_RegionH = ((sin(pi.*r).*(Vin + Vo.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) - Vin.*cos(pi.*Dy1.*r) - Vo.*cos((r.*pi.*(Dy1 - 2.*Dp + Dy2 + 2))/2)))/Z0 + ((cos(pi.*r) + 1).*(Vo.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) + Vin.*sin(pi.*Dy1.*r) + Vo.*sin((r.*pi.*(Dy1 - 2.*Dp + Dy2 + 2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));
        I_zvs_s34_RegionH = -(((cos(pi.*r) + 1).*(Vin.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2))/2) - Vin.*sin(pi.*r) + Vo.*sin(pi.*Dy2.*r) + Vin.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2)))/Z0 + (sin(pi.*r).*(Vin + Vo - Vin.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2))/2) + Vin.*cos(pi.*r) - Vo.*cos(pi.*Dy2.*r) - Vin.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));

        Imax2 = abs((Vin.^2.*cos(pi.*r.*(Dy1 - 1)) + Vo.^2.*cos(pi.*r.*(Dy2 - 1)) + Vin.^2 + Vo.^2 - Vin.*Vo.*cos((pi.*r.*(2.*Dp + Dy1 + Dy2 - 2))/2) - Vin.*Vo.*cos((pi.*r.*(2.*Dp + Dy1 - Dy2))/2) - Vin.*Vo.*cos((pi.*r.*(2.*Dp - Dy1 + Dy2))/2) - Vin.*Vo.*cos((pi.*r.*(2.*Dp - Dy1 - Dy2 + 2))/2))./(Z0.^2.*(cos(pi.*r) + 1)));
        Imax1 = abs((2.*Vin.^2 - 2.*Vo.^2.*cos((Dy2.*pi.*r)/2).^2 - 2.*Vin.^2.*cos((Dy1.*pi.*r)/2).^2 + 2.*Vo.^2 + 4.*Vin.*Vo.*cos(pi.*r).*cos(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2) + 4.*Vin.*Vo.*sin(pi.*r).*sin(Dp.*pi.*r).*sin((Dy1.*pi.*r)/2).*sin((Dy2.*pi.*r)/2))./(Z0.^2.*(cos(pi.*r) + 1)));

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


function selected_vars = select_pareto_solution(optimal_vars, Irms_values, PoN_values, DoF)
    % SELECT_PARETO_SOLUTION - Interactive tool to select a specific point from the Pareto front
    %
    % Inputs:
    %   optimal_vars - Matrix of optimization variables for each Pareto point
    %   Irms_values - Vector of Irms values for each Pareto point
    %   PoN_values - Vector of PoN values for each Pareto point
    %   DoF - Degree of freedom (number of optimization variables)
    %
    % Outputs:
    %   selected_vars - The selected optimization variables
    
        % Plot Pareto front
        fig = figure('Name', 'Select Pareto Solution', 'Position', [100, 100, 800, 600]);
        scatter(PoN_values, Irms_values, 5, 'filled');
        hold on;
        
        % Add labels
        xlabel('Normalized output power PoN', 'Interpreter', 'latex');
        ylabel('Normalized tank current I_{rms}', 'Interpreter', 'latex');
        title('Click to select a solution from the Pareto front', 'Interpreter', 'latex');
        grid on;
        
        % Instructions
        disp('Click on a point in the Pareto front to select that solution.');
        
        % Wait for user to click on a point
        [x, y] = ginput(1);
        
        % Find the closest point on the Pareto front
        distances = sqrt((PoN_values - x).^2 + (Irms_values - y).^2);
        [~, idx] = min(distances);
        
        % Highlight the selected point
        scatter(PoN_values(idx), Irms_values(idx), 100, 'r', 'filled');
        
        % Display the selected solution
        selected_vars = optimal_vars(idx,:);
        
        % Create text annotation with the details based on DoF
        switch DoF
            case 1
                annotation('textbox', [0.15, 0.15, 0.3, 0.2], 'String', ...
                    {sprintf('Selected solution:'), ...
                     sprintf('PoN = %.4f', PoN_values(idx)), ...
                     sprintf('Irms = %.4f', Irms_values(idx)), ...
                     sprintf('Dp = %.4f', selected_vars(1))}, ...
                    'FitBoxToText', 'on', 'BackgroundColor', 'white');
                    
                disp(['Selected Point: PoN = ', num2str(PoN_values(idx)), ', Irms = ', num2str(Irms_values(idx))]);
                disp(['Variables: Dp = ', num2str(selected_vars(1))]);
                
            case 2
                annotation('textbox', [0.15, 0.15, 0.3, 0.2], 'String', ...
                    {sprintf('Selected solution:'), ...
                     sprintf('PoN = %.4f', PoN_values(idx)), ...
                     sprintf('Irms = %.4f', Irms_values(idx)), ...
                     sprintf('Dp = %.4f', selected_vars(1)), ...
                     sprintf('Dy1/Dy2 = %.4f', selected_vars(2))}, ...
                    'FitBoxToText', 'on', 'BackgroundColor', 'white');
                    
                disp(['Selected Point: PoN = ', num2str(PoN_values(idx)), ', Irms = ', num2str(Irms_values(idx))]);
                disp(['Variables: Dp = ', num2str(selected_vars(1)), ...
                      ', Dy1/Dy2 = ', num2str(selected_vars(2))]);
                      
            case 3
                annotation('textbox', [0.15, 0.15, 0.3, 0.2], 'String', ...
                    {sprintf('Selected solution:'), ...
                     sprintf('PoN = %.4f', PoN_values(idx)), ...
                     sprintf('Irms = %.4f', Irms_values(idx)), ...
                     sprintf('Dp = %.4f', selected_vars(1)), ...
                     sprintf('Dy1 = %.4f', selected_vars(2)), ...
                     sprintf('Dy2 = %.4f', selected_vars(3))}, ...
                    'FitBoxToText', 'on', 'BackgroundColor', 'white');
                    
                disp(['Selected Point: PoN = ', num2str(PoN_values(idx)), ', Irms = ', num2str(Irms_values(idx))]);
                disp(['Variables: Dp = ', num2str(selected_vars(1)), ...
                      ', Dy1 = ', num2str(selected_vars(2)), ...
                      ', Dy2 = ', num2str(selected_vars(3))]);
                      
            case 4
                annotation('textbox', [0.15, 0.15, 0.3, 0.2], 'String', ...
                    {sprintf('Selected solution:'), ...
                     sprintf('PoN = %.4f', PoN_values(idx)), ...
                     sprintf('Irms = %.4f', Irms_values(idx)), ...
                     sprintf('Dp = %.4f', selected_vars(1)), ...
                     sprintf('Dy1 = %.4f', selected_vars(2)), ...
                     sprintf('Dy2 = %.4f', selected_vars(3)), ...
                     sprintf('r = %.4f', selected_vars(4))}, ...
                    'FitBoxToText', 'on', 'BackgroundColor', 'white');
                    
                disp(['Selected Point: PoN = ', num2str(PoN_values(idx)), ', Irms = ', num2str(Irms_values(idx))]);
                disp(['Variables: Dp = ', num2str(selected_vars(1)), ...
                      ', Dy1 = ', num2str(selected_vars(2)), ...
                      ', Dy2 = ', num2str(selected_vars(3)), ...
                      ', r = ', num2str(selected_vars(4))]);
        end
    end