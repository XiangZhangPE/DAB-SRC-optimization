% Constants
fs = 290e3; % Switching frequency in Hz
Lr = 12e-6; % Resonant inductor in H
Cr = 48e-9; % Resonant capacitor in F
Z0 = sqrt(Lr/Cr); % Characteristic impedance
fr = 1/(2*pi*sqrt(Lr*Cr)); % Resonant frequency in Hz
td = 50e-9; % Dead time in s
Dd = td*fs; % Duty cycle of dead time
r = fr/fs; % Ratio of resonant frequency to switching frequency
omega_r = 2*pi*fr; % Angular resonant frequency
Ln = 100; % Inductance normalization factor

Vin = 800; % Input voltage in V
Vo = 1000; % Output voltage in V

ILm0 = pi*r*Vo/(2*Ln*Z0); % Initial magnetizing current
Ioss = 0.06; % Output current

% Create a meshgrid for duty cycles
[Dy1, Dy2, Dp] = meshgrid(linspace(0, 1, 100), linspace(0, 1, 100), linspace(0, 1, 100));
hand = figure; % Create a new figure

% Calculate switching currents in different stages
IL0_RegionE = -(((cos(pi.*r) + 1).*(Vo.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) + Vo.*sin(pi.*r) + Vin.*sin(pi.*Dy1.*r) - Vo.*sin((r.*pi.*(Dy1 - 2.*Dp + Dy2))/2)))/Z0 + (sin(pi.*r).*(Vin - Vo + Vo.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2 - 2))/2) - Vo.*cos(pi.*r) - Vin.*cos(pi.*Dy1.*r) + Vo.*cos((r.*pi.*(Dy1 - 2.*Dp + Dy2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));
IL1_RegionE = -(Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2) + Vin.*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2) + Vin.*sin(pi.*r.*(Dy1 - 1)) + Vo.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2).*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 4))/2).*cos(pi.*r) + Vin.*cos(pi.*r.*(Dy1 - 1)).*sin(pi.*r) + Vin.*sin(pi.*r.*(Dy1 - 1)).*cos(pi.*r) + Vo.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2).*sin(pi.*r) + Vo.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2 - 2))/2).*cos(pi.*r))./(2.*Z0.*(cos(pi.*r) + 1));
IL2_RegionE = (((cos(pi.*r) + 1).*(Vin.*sin((r.*pi.*(2.*Dp - Dy1 + Dy2))/2) - Vin.*sin(pi.*r) + Vo.*sin(pi.*Dy2.*r) + Vin.*sin((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2)))/Z0 + (sin(pi.*r).*(Vin + Vo - Vin.*cos((r.*pi.*(2.*Dp - Dy1 + Dy2))/2) + Vin.*cos(pi.*r) - Vo.*cos(pi.*Dy2.*r) - Vin.*cos((r.*pi.*(2.*Dp + Dy1 + Dy2 - 2))/2)))/Z0)./(2.*(cos(pi.*r) + 1));
IL3_RegionE = (Vin.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2))/2) + Vo.*sin(pi.*r) + Vin.*sin((r.*pi.*(2.*Dp - Dy1 - Dy2 + 2))/2) + Vo.*sin(pi.*r.*(Dy2 - 1)) + Vin.*cos(pi.*r).*sin((r.*pi.*(2.*Dp - Dy1 - Dy2 + 2))/2) - Vin.*sin(pi.*r).*cos((r.*pi.*(2.*Dp - Dy1 - Dy2 + 2))/2) + Vo.*cos(pi.*r.*(Dy2 - 1)).*sin(pi.*r) + Vo.*sin(pi.*r.*(Dy2 - 1)).*cos(pi.*r) - Vin.*cos((r.*pi.*(2.*Dp + Dy1 - Dy2))/2).*sin(pi.*r) + Vin.*sin((r.*pi.*(2.*Dp + Dy1 - Dy2))/2).*cos(pi.*r))./(2.*Z0.*(cos(pi.*r) + 1));

% Define zero-voltage switching (ZVS) conditions
I_zvs_p12_RegionE =  IL1_RegionE <=0;
I_zvs_s12_RegionE =  IL3_RegionE >=0;
I_zvs_p34_RegionE =  -IL0_RegionE >=0;
I_zvs_s34_RegionE =  -IL2_RegionE <=0;

% Optimization problem setup
N = 100;
Dy1 = optimvar('Dy1',N,'LowerBound',0,'UpperBound',1);
Dy2 = optimvar('Dy2',N,'LowerBound',0,'UpperBound',1);
Dp = optimvar('Dp',N,'LowerBound',0,'UpperBound',0.5);
minIrms_prob = optimproblem;

% Define constraints for the optimization problem
minIrms_prob.Constraints.cons1 = I_zvs_p12_RegionE;
minIrms_prob.Constraints.cons2 = I_zvs_s12_RegionE;
minIrms_prob.Constraints.cons3 = I_zvs_p34_RegionE;
minIrms_prob.Constraints.cons4 = I_zvs_s34_RegionE;

% Define the objective function for the optimization problem
minIrms_prob.Objective = Irms_func; % Set the objective function

% Initial guess for the optimization variables
init.Dy1 = 0.5;
init.Dy2 = 0.5;
init.Dp = 0.3;

init.r = r;
init.Vo = Vo;
init.Vin = Vin;
init.Z0 = Z0;

% Solve the optimization problem
[sol,fval,exitflag,output] = solve(minIrms_prob,init);

% Extract the optimized values
Dy1_opt = sol.Dy1;
Dy2_opt = sol.Dy2;
Dp_opt = sol.Dp;

% Plot the solution
figure(hand)
plot3(Dy1_opt, Dy2_opt, Dp_opt, 'r.', 'MarkerSize', 25)
xlabel('Dy1')
ylabel('Dy2')
zlabel('Dp')
title('Optimization Results')
grid on
hold off









% Calculate power normalized function
function f = PoN_func(Dy1,Dy2,Dp,r,Vo,Vin,Z0)
f = ((cos(pi.*r.*(Dp-0.5)).*cos(pi.*r.*(1-Dy1)./2).*cos(pi.*r.*(1-Dy2)./2))./(2.*r.*cos(pi.*r./2))-1./(2.*r));
end
% Calculate Irms normalized function
function g = Irms_func(Dy1,Dy2,Dp,r,Vo,Vin,Z0)
g = Z0./Vin.*sqrt((1/2/pi./r).*(-(2.*Vin.^2.*sin(Dy1.*pi.*r) + 2.*Vo.^2.*sin(Dy2.*pi.*r) + 2.*Vin.^2.*sin(pi.*r.*(Dy1 - 1)) + 2.*Vo.^2.*sin(pi.*r.*(Dy2 - 1)) + 2.*Vin.^2.*sin(pi.*r) + 2.*Vo.^2.*sin(pi.*r) + 2.*Vin.*Vo.*sin((pi.*r.*(2.*Dp - Dy1 - Dy2 + 2))./2) - 2.*Vin.*Vo.*sin((pi.*r.*(Dy1 - 2.*Dp + Dy2))./2) + 2.*Vin.*Vo.*sin((pi.*r.*(2.*Dp + Dy1 + Dy2 - 2))./2) + 2.*Vin.*Vo.*sin((pi.*r.*(2.*Dp + Dy1 + Dy2 - 4))./2) + 2.*Vin.*Vo.*sin((pi.*r.*(2.*Dp + Dy1 - Dy2))./2) + 2.*Vin.*Vo.*sin((pi.*r.*(2.*Dp - Dy1 + Dy2))./2) + 2.*Vin.*Vo.*sin((pi.*r.*(2.*Dp + Dy1 - Dy2 - 2))./2) + 2.*Vin.*Vo.*sin((pi.*r.*(2.*Dp - Dy1 + Dy2 - 2))./2) - 2.*Vin.^2.*pi.*r - 2.*Vo.^2.*pi.*r + 2.*Vin.^2.*pi.*r.*cos(Dy1.*pi.*r) + 2.*Vo.^2.*pi.*r.*cos(Dy2.*pi.*r) - 2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(Dy1 - 2.*Dp + Dy2))./2) + 4.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp + Dy1 + Dy2 - 2))./2) + 2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp + Dy1 + Dy2 - 4))./2) - 2.*Dy1.*Vin.^2.*pi.*r.*cos(Dy1.*pi.*r) - 2.*Dy2.*Vo.^2.*pi.*r.*cos(Dy2.*pi.*r) + 2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp + Dy1 - Dy2))./2) + 2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp - Dy1 + Dy2))./2) - 2.*Dy1.*Vin.^2.*pi.*r.*cos(pi.*r.*(Dy1 - 1)) - 2.*Dy2.*Vo.^2.*pi.*r.*cos(pi.*r.*(Dy2 - 1)) - 2.*Dp.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp - Dy1 - Dy2 + 2))./2) + Dy1.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp - Dy1 - Dy2 + 2))./2) + Dy2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp - Dy1 - Dy2 + 2))./2) - 2.*Dp.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(Dy1 - 2.*Dp + Dy2))./2) - 2.*Dp.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp + Dy1 + Dy2 - 2))./2) - 2.*Dp.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp + Dy1 + Dy2 - 4))./2) + Dy1.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(Dy1 - 2.*Dp + Dy2))./2) + Dy2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(Dy1 - 2.*Dp + Dy2))./2) - Dy1.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp + Dy1 + Dy2 - 2))./2) - Dy2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp + Dy1 + Dy2 - 2))./2) - Dy1.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp + Dy1 + Dy2 - 4))./2) - Dy2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp + Dy1 + Dy2 - 4))./2) - 2.*Dp.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp + Dy1 - Dy2))./2) - 2.*Dp.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp - Dy1 + Dy2))./2) - 2.*Dp.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp + Dy1 - Dy2 - 2))./2) - 2.*Dp.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp - Dy1 + Dy2 - 2))./2) - Dy1.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp + Dy1 - Dy2))./2) + Dy1.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp - Dy1 + Dy2))./2) + Dy2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp + Dy1 - Dy2))./2) - Dy2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp - Dy1 + Dy2))./2) - Dy1.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp + Dy1 - Dy2 - 2))./2) + Dy1.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp - Dy1 + Dy2 - 2))./2) + Dy2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp + Dy1 - Dy2 - 2))./2) - Dy2.*Vin.*Vo.*pi.*r.*cos((pi.*r.*(2.*Dp - Dy1 + Dy2 - 2))./2))./(4.*Z0.^2.*(cos(pi.*r) + 1))));
end
