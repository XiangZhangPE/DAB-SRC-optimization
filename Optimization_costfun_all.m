clc, close all;
%% Define constants
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

% plot background ZVS figure
handle = figure(3);
plotoptions = 2;
% r = 1/1.8;
% ZVSPlot_DABsrc(r,Vo,Vin,Z0, plotoptions); hold on;
% r = 1/1.5;
% ZVSPlot_DABsrc(r,Vo,Vin,Z0, plotoptions); hold on;
% r = 1/1.1;
% ZVSPlot_DABsrc(r,Vo,Vin,Z0, plotoptions); hold on;
% r = fs/fr;
ZVSPlot_DABsrc(r,Vo,Vin,Z0, plotoptions); hold on;

% Define optimization problem options
options = optimoptions('fmincon', 'Display', 'off');
% options = optimoptions('fmincon', ...
%     'OptimalityTolerance', 1e-6, ...  % improve convergence accurancy
%     'FunctionTolerance', 1e-6, ...
%     'StepTolerance', 1e-10, ...
%     'ConstraintTolerance', 1e-9, ...
%     'Display', 'off');  % 观察每次迭代结果

%% low high power guess for mini Irms

% bounds = [0 0; 0.5 1];
% DoF = 2; % Degree of freedom
% initial_guess = [0.06, 0.7];  % Initial guess
% % PoN_values = linspace(0.015, 0.45, 100);
% PoN_values = linspace(0.015, 0.3, 100);

% bounds = [0 0 0; 0.5 1 1]; % Lower and upper bounds
% DoF = 3; % Degree of freedom
% initial_guess = [0.06, 0.7, 0.95];  % Initial guess
% PoN_values = linspace(0.015, 0.45, 100);

% % M<1 (Ioff lambad=400) (Irms lambad=500)
% bounds = [0 0 0 1/1.2; 0.5 1 1 1/1.05]; % Lower and upper bounds
% DoF = 4; % Degree of freedom
% initial_guess = [0.06, 0.7, 0.95, 1/1.1];  % Initial guess
% % PoN_values = linspace(0.015, 0.3, 100);
% PoN_values = linspace(0.33, 0.45, 100);

% % M>1 (Ioff lambad=400) (Irms lambad=500)
% bounds = [0 0 0 1/1.2; 0.5 1 1 1/1.1]; % Lower and upper bounds
% DoF = 4; % Degree of freedom
% initial_guess = [0.06, 0.95, 0.7, 1/1.1];  % Initial guess
% PoN_values = linspace(0.014, 0.3, 100);

% % M=580/800, fs = 250e3 Irms lambad=1000
% bounds = [0 0 0; 0.5 1 1];
% DoF = 3; % Degree of freedom
% initial_guess = [0.06, 0.7, 0.95];  % Initial guess
% PoN_values = linspace(0.015, 0.45, 100);

%% meduim to high power guess for mini Irms
% % TPS optimize ABE
% bounds = [0 0 0; 0.5 1 1]; % Lower and upper bounds
% DoF = 3; % Degree of freedom
% initial_guess = [0.25, 1, 1];  % Initial guess
% PoN_values = linspace(0.45, 1.2, 100);

% % TPS optimize EFGH
% bounds = [0 0 0; 1 1 1]; % Lower and upper bounds
% DoF = 3; % Degree of freedom
% initial_guess = [0.35, 1, 1];  % Initial guess
% PoN_values = linspace(0.05, 1.2, 100);

% % VFM+TPS optimize ABE lambda = 400
% bounds = [0 0 0 1/1.2; 0.5 1 1 1/1.01]; % Lower and upper bounds
% DoF = 4; % Degree of freedom
% initial_guess = [0.25, 1, 1, 1/1.2];  % Initial guess
% % initial_guess = [0.35, 1, 1, 1/1.15];  % Initial guess
% PoN_values = linspace(0.45, 1.2, 100);

%% meduim to high power guess for mini Ioff
% bounds = [0 0 0; 0.5 1 1]; % Lower and upper bounds
% DoF = 3; % Degree of freedom
% initial_guess = [0.25, 1, 1];  % Initial guess
% PoN_values = linspace(0.45, 1.2, 100);

bounds = [0 0 0 1/1.2; 0.5 1 1 1/1.1]; % Lower and upper bounds
DoF = 4; % Degree of freedom
initial_guess = [0.5, 1, 1, 1/1.2];  % Initial guess
PoN_values = linspace(0.45, 1.2, 100);

%% all power guess
% bounds = [0 0 0 1/1.2; 0.5 1 1 1]; % Lower and upper bounds
% DoF = 4; % Degree of freedom
% initial_guess = [0.06, 0.7, 0.95, 1/1.1];  % Initial guess
% PoN_values = linspace(0.02, 1.1, 100);

%% Run optimization for each PoN target value

optimal_points = zeros(length(PoN_values), DoF);
optimal_power = zeros(length(PoN_values), 1);
optimal_Irms = zeros(length(PoN_values), 1);
Ioffp12 = zeros(length(PoN_values), 1);
Ioffp34 = zeros(length(PoN_values), 1);
Ioffs12 = zeros(length(PoN_values), 1);
Ioffs34 = zeros(length(PoN_values), 1);

for i = 1:length(PoN_values)
    PoN_target = PoN_values(i);
    
    % Define objective function with penalty
    obj_fun = @(vars) J_function(DoF, vars, Vin, Vo, PoN_target);  % objective function now include PoN_target
    
    % Define nonlinear constraints, import iteration vars, Vin, Vo, Z0
    nonlcon = @(vars) zvs_function(DoF, vars, Vin, Vo, Z0);
    
    % Run optimization
    [optimal_vars, fval, exitflag, output] = fmincon(obj_fun, initial_guess, [], [], [], [], bounds(1,:), bounds(2,:), nonlcon, options);

    % Store result
    optimal_points(i,:) = optimal_vars;
    optimal_power(i,:) = f_function(DoF, optimal_vars, Vin, Vo);
    optimal_Irms(i,:) = g_function(DoF, optimal_vars, Vin, Vo);
    [Ip12, Ip34, Is12, Is34, ILmax] = h_function(DoF, optimal_vars, Vin, Vo);
    Ioffp12(i,:) = Ip12;
    Ioffp34(i,:) = Ip34;
    Ioffs12(i,:) = Is12;
    Ioffs34(i,:) = Is34;
    ILmax(i,:) = ILmax;
    
    % Update initial guess for next optimization
    initial_guess = optimal_vars;
    
    % output virsialization and optimization status
    if exitflag > 0
        fprintf('✅ Optimization successful for PoN_target = %.4f. Final Cost: %.6f\n', PoN_target, fval);
        if DoF >= 4
            fprintf('Optimal Variables: Dp = %.6f, Dy1 = %.6f, Dy2 = %.6f\n, r = %.6f\n', optimal_vars(1), optimal_vars(2), optimal_vars(3), optimal_vars(4));
        elseif DoF == 3
            fprintf('Optimal Variables: Dp = %.6f, Dy1 = %.6f, Dy2 = %.6f\n', optimal_vars(1), optimal_vars(2), optimal_vars(3));
        elseif DoF == 2
            fprintf('Optimal Variables: Dp = %.6f, Dy1(Dy2) = %.6f\n', optimal_vars(1), optimal_vars(2));
        else
            fprintf('Optimal Variable: Dp = %.6f\n', optimal_vars(1));
        end
    else
        warning('⚠️ Optimization did NOT converge for PoN_target = %.4f! Exit flag: %d\n', PoN_target, exitflag);
        fprintf('Reason: %s\n', output.message);
        if DoF >= 4
            fprintf('Optimal Variables: Dp = %.6f, Dy1 = %.6f, Dy2 = %.6f\n, r = %.6f\n', optimal_vars(1), optimal_vars(2), optimal_vars(3), optimal_vars(4));
        elseif DoF == 3
            fprintf('Last Attempt Variables: Dp = %.6f, Dy1 = %.6f, Dy2 = %.6f\n', optimal_vars(1), optimal_vars(2), optimal_vars(3));
        elseif DoF == 2
            fprintf('Last Attempt Variables: Dp = %.6f, Dy1(Dy2) = %.6f\n', optimal_vars(1), optimal_vars(2));
        else
            fprintf('Last Attempt Variable: Dp = %.6f\n', optimal_vars(1));
        end
    end
end

figure(1)
% plot(PoN_values, optimal_Irms, 'b-', 'LineWidth', 1); hold on;
% ylabel('Normalized tank current Irms', 'Interpreter', 'latex');
% title('Irms under different Power','Interpreter', 'latex');

plot(PoN_values, Ioffp12, 'r-', 'LineWidth', 1, 'DisplayName', 'I_{off,P12}'); hold on;
plot(PoN_values, Ioffp34, 'b-', 'LineWidth', 1, 'DisplayName', 'I_{off,P34}'); hold on;
plot(PoN_values, Ioffs12, 'm-', 'LineWidth', 1, 'DisplayName', 'I_{off,S12}'); hold on;
plot(PoN_values, Ioffs34, 'c-', 'LineWidth', 1, 'DisplayName', 'I_{off,S34}'); hold on;
plot(PoN_values, ILmax, 'k--', 'LineWidth', 1, 'DisplayName', 'I_{L,max}'); hold on;
ylabel('Normalized turn-off current', 'Interpreter', 'latex');
title('Ioff under different Power','Interpreter', 'latex');

xlabel('Normalized output power PoN', 'Interpreter', 'latex');

legend('show');
% xlim([0 0.45]);
% ylim([0 0.5]);
grid on;


figure(2)
plot(PoN_values, optimal_points(:,1), 'k-', 'LineWidth', 1, 'DisplayName', 'Dp'); hold on;
plot(PoN_values, optimal_points(:,2), 'r-', 'LineWidth', 1, 'DisplayName', 'Dy1'); hold on;
plot(PoN_values, optimal_points(:,3), 'b-', 'LineWidth', 1, 'DisplayName', 'Dy2'); hold on;
xlabel('Normalized output power PoN', 'Interpreter', 'latex');
ylabel('Phase-shift ratios','Interpreter', 'latex');
legend('show');
% xlim([0 0.45]);
% ylim([0 0.5]);
grid on;

figure(3)
if DoF >= 3
    % Visualize results in 3D
    scatter3(optimal_points(:,2), optimal_points(:,3), optimal_points(:,1), 10, 'red','filled'); % 50 is the point size
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





% Define Cost function with penalty for PoN_target
function J = J_function(DoF, vars, Vin, Vo, PoN_target)

    % calculate PoN
    PoN = f_function(DoF, vars, Vin, Vo);
    % calculate Irms
    Irms = g_function(DoF, vars, Vin, Vo);
    Ioff = h_function(DoF, vars, Vin, Vo);
    
    % construct cost function
    lambda = 500;  % add panelty coefficient (weight λ adjustable)
    penalty = lambda * (PoN - PoN_target)^2;
    
    % final cost function
    J = Irms + penalty;
end

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