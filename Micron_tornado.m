%% Sensitivity Analysis: Differences from Baseline (Tornado Plot)
% Parameters:
% Temperature: Baseline = 873 K, Varied by ±200 K
% Pressure: Baseline = 66.6612 Pa, Varied by ±40%
% H2 Concentration: Baseline = 10 mol/m^3, Varied by ±40%
% Delta: Baseline = 0.5 m, Varied by ±40%
% Gas Inlet (SiH4 Concentration): Baseline = 1000 mol/m^3, Varied by ±40%

%% Define Baseline Values
T_base    = 273 + 600;       % 873 K
P_base    = 66.6612;         % Pa
H2_base   = 10;              % mol/m^3
delta_base = 500e-3;         % 0.5 m
SiH4_base = 1000;            % mol/m^3

%% Define Variation Ranges
% Temperature variation: ±200 K
T_low  = 580+273;       % 673 K
T_high = 725+273;       % 1073 K

% Pressure variation: ±40%
variation = 0.40;
P_low  = 1;
P_high = 100;

% H2 concentration variation: ±40%
H2_low  = 100;
H2_high = 1000;

% Delta variation: ±40%
delta_low  = delta_base * (1 - variation);
delta_high = delta_base * (1 + variation);

% Gas inlet (SiH4 concentration) variation: ±40%
SiH4_low  = SiH4_base * (1 - variation);
SiH4_high = 1000;

%% Run Baseline Simulation (Extended Time Span [0,2000] s)
[t_base, thickness_base] = runSimulation(T_base, P_base, H2_base, delta_base, SiH4_base);
time_base = findTimeToReachThickness(t_base, thickness_base, 1000);

%% Run Simulations for Each Parameter Variation and Compute Differences from Baseline

% Temperature Variations (other parameters at baseline)
[t_T_low, thickness_T_low] = runSimulation(T_low, P_base, H2_base, delta_base, SiH4_base);
time_T_low = findTimeToReachThickness(t_T_low, thickness_T_low, 1000);
[t_T_high, thickness_T_high] = runSimulation(T_high, P_base, H2_base, delta_base, SiH4_base);
time_T_high = findTimeToReachThickness(t_T_high, thickness_T_high, 1000);
dev_T_low  = time_T_low  - time_base;
dev_T_high = time_T_high - time_base;

% Pressure Variations
[t_P_low, thickness_P_low] = runSimulation(T_base, P_low, H2_base, delta_base, SiH4_base);
time_P_low = findTimeToReachThickness(t_P_low, thickness_P_low, 1000);
[t_P_high, thickness_P_high] = runSimulation(T_base, P_high, H2_base, delta_base, SiH4_base);
time_P_high = findTimeToReachThickness(t_P_high, thickness_P_high, 1000);
dev_P_low  = time_P_low  - time_base;
dev_P_high = time_P_high - time_base;

% H2 Concentration Variations
[t_H2_low, thickness_H2_low] = runSimulation(T_base, P_base, H2_low, delta_base, SiH4_base);
time_H2_low = findTimeToReachThickness(t_H2_low, thickness_H2_low, 1000);
[t_H2_high, thickness_H2_high] = runSimulation(T_base, P_base, H2_high, delta_base, SiH4_base);
time_H2_high = findTimeToReachThickness(t_H2_high, thickness_H2_high, 1000);
dev_H2_low  = time_H2_low  - time_base;
dev_H2_high = time_H2_high - time_base;

% Delta Variations
[t_delta_low, thickness_delta_low] = runSimulation(T_base, P_base, H2_base, delta_low, SiH4_base);
time_delta_low = findTimeToReachThickness(t_delta_low, thickness_delta_low, 1000);
[t_delta_high, thickness_delta_high] = runSimulation(T_base, P_base, H2_base, delta_high, SiH4_base);
time_delta_high = findTimeToReachThickness(t_delta_high, thickness_delta_high, 1000);
dev_delta_low  = time_delta_low  - time_base;
dev_delta_high = time_delta_high - time_base;

% Gas Inlet (SiH4 Concentration) Variations
[t_SiH4_low, thickness_SiH4_low] = runSimulation(T_base, P_base, H2_base, delta_base, SiH4_low);
time_SiH4_low = findTimeToReachThickness(t_SiH4_low, thickness_SiH4_low, 1000);
[t_SiH4_high, thickness_SiH4_high] = runSimulation(T_base, P_base, H2_base, delta_base, SiH4_high);
time_SiH4_high = findTimeToReachThickness(t_SiH4_high, thickness_SiH4_high, 1000);
dev_SiH4_low  = time_SiH4_low  - time_base;
dev_SiH4_high = time_SiH4_high - time_base;

%% Prepare Data for the Tornado Plot (Differences from Baseline)
params    = {'Temperature', 'Pressure', 'H2 Concentration', 'Delta', 'SiH4 Inlet'};
% For each parameter, the "low" and "high" differences from baseline:
devLow  = [dev_T_low,  dev_P_low,  dev_H2_low,  dev_delta_low,  dev_SiH4_low];
devHigh = [dev_T_high, dev_P_high, dev_H2_high, dev_delta_high, dev_SiH4_high];

%% Create the Tornado Plot (Centered at 0 = Baseline)
figure;
hold on;
for i = 1:length(params)
    % Plot error bars centered at 0.
    % The left error is -devLow (if devLow is negative) and the right error is devHigh.
    errorbar(0, i, -devLow(i), devHigh(i), 'horizontal', 'o', 'LineWidth', 2);
end
set(gca, 'YTick', 1:length(params), 'YTickLabel', params);
xlabel('Deviation from Baseline Time (s)');
title('Tornado Plot: Deviations from Baseline');
grid on;
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper Functions

function [t, thickness_A] = runSimulation(T_mod, P_mod, H2_mod, delta_mod, SiH4_mod)
    % Load baseline constants from Micron_constants1.m
    [R, ~, E_decom, A_decom, SiH4_conc, ~, tspan, ~, delta, ...
     D_SiH4, D_SiH2, D_H2, A_SiH4, E_SiH4, A_SiH2, E_SiH2, A_H2des, E_H2des, ...
     A_SiH2_growth, E_SiH2_growth, Rad, A, Thick, DS, MMs, Ratio, A_tot, ...
     A_dep_SiH4, E_dep_SiH4, A_dep_growth, E_dep_growth] = Micron_constants1();
    
    % Override baseline values with modified inputs:
    T = T_mod;
    P = P_mod;
    H2_conc = H2_mod;
    delta = delta_mod;
    SiH4_conc = SiH4_mod;
    
    % Use an extended time span to capture slower growth scenarios
    tspan = [0, 2000];
    
    % Recalculate diffusivities (which depend on T)
    D_SiH4 = (1/100^2) * (9.6e-5) * (T^1.5);
    D_SiH2 = (1/100^2) * (3.45e-5) * (T^1.5);
    D_H2   = (1/100^2) * (3.45e-5) * (T^1.5);
    
    % Initial conditions: [SiH4_g; SiH2_g; H2_g; SiH4_b; SiH2_b; H2_b]
    C0 = [SiH4_conc; 0; H2_conc; 0; 0; 0];
    
    % Reaction constant via the Arrhenius relation
    k1 = A_decom * exp(-E_decom/(R*T)) * 0.01;
    
    % Solve the ODE system
    [t, C] = ode45(@(t, C) silaneDiffusionODE(t, C, k1, R, T, P, delta, D_SiH4, D_SiH2, D_H2), tspan, C0);
    
    % Calculate deposition flux and film growth rate
    k_SiH4_dep = A_dep_SiH4 * exp(-E_dep_SiH4/(R*T));
    k_p_dep    = A_dep_growth * exp(-E_dep_growth/(R*T));
    J = (k_p_dep * k_SiH4_dep .* C(:,4)) ./ (1 + k_SiH4_dep .* C(:,4));
    growthRate_m_s = (J * 28.085) / 2330;
    growthRate_A_s = growthRate_m_s * 1e10;  % Convert m/s to Å/s
    thickness_A = cumtrapz(t, growthRate_A_s); % Integrated film thickness in Å
    
    % Nested ODE function
    function dCdt = silaneDiffusionODE(~, C, k1, R, T, P, delta, D_SiH4, D_SiH2, D_H2)
        SiH4_g = C(1);
        SiH2_g = C(2);
        H2_g   = C(3);
        SiH4_b = C(4);
        SiH2_b = C(5);
        H2_b   = C(6);
        
        dCdt = zeros(6,1);
        % Gas-phase silane decomposition
        dCdt(1) = -k1 * SiH4_g;
        dCdt(2) =  k1 * SiH4_g;
        dCdt(3) =  2 * k1 * SiH4_g;
        % Diffusion (Fick's law)
        dCdt(4) = (D_SiH4/delta) * (SiH4_g - SiH4_b);
        dCdt(5) = (D_SiH2/delta) * (SiH2_g - SiH2_b);
        dCdt(6) = (D_H2/delta)   * (H2_g   - H2_b);
    end
end

function time_to_thickness = findTimeToReachThickness(t, thickness_A, targetThickness)
    idx = find(thickness_A >= targetThickness, 1, 'first');
    if ~isempty(idx)
        time_to_thickness = t(idx);
    else
        time_to_thickness = NaN;
    end
end
