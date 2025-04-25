%% Define Baseline Values
T_base     = 273 + 600;       % 873 K
P_base     = 66.6612;         % Pa
H2_base    = 10;              % SCCM
delta_base = 500e-3;          % 0.5 m
SiH4_base  = 1000;            % SCCM

%% Define Variation Ranges
% Temperature variation: ±200 K
T_low  = 580;       % 673 K
T_high = 725;       % 1073 K

% Pressure variation: ±40%
variation = 0.40;
P_low  = 1;
P_high = 100;

% H2 concentration variation: ±40%
H2_low  = 100;
H2_high = 100;

% Delta variation: ±40%
delta_low  = .001;
delta_high = .005;

% Gas inlet (SiH4 concentration) variation: ±40%
SiH4_low  = SiH4_base * (1 - variation);
SiH4_high = 1000;

%% Run Baseline Simulation
[t_base, thickness_base] = runSimulation(T_base, P_base, H2_base, delta_base, SiH4_base);
time_base = findTimeToReachThickness(t_base, thickness_base, 1000);

%% Grid Search Over All Parameter Combinations
% Candidate arrays (each has two entries: low and high)
T_candidates     = [T_low, T_high];
P_candidates     = [P_low, P_high];
H2_candidates    = [H2_low, H2_high];
delta_candidates = [delta_low, delta_high];
SiH4_candidates  = [SiH4_low, SiH4_high];

nComb = 2^5;  % Total number of combinations (32)
results = zeros(nComb, 6);  % Columns: T, P, H2, delta, SiH4, time

combIndex = 1;
for iT = 1:length(T_candidates)
    for iP = 1:length(P_candidates)
        for iH2 = 1:length(H2_candidates)
            for idelta = 1:length(delta_candidates)
                for iSiH4 = 1:length(SiH4_candidates)
                    % Get current combination
                    T_val     = T_candidates(iT);
                    P_val     = P_candidates(iP);
                    H2_val    = H2_candidates(iH2);
                    delta_val = delta_candidates(idelta);
                    SiH4_val  = SiH4_candidates(iSiH4);
                    
                    % Run simulation with extended time span [0,2000] s
                    [t_sim, thickness_sim] = runSimulation(T_val, P_val, H2_val, delta_val, SiH4_val);
                    time_val = findTimeToReachThickness(t_sim, thickness_sim, 1000);
                    
                    results(combIndex,:) = [T_val, P_val, H2_val, delta_val, SiH4_val, time_val];
                    combIndex = combIndex + 1;
                end
            end
        end
    end
end

%% Find the Combination with the Minimum Time (ignoring NaNs)
valid_results = results(~isnan(results(:,6)), :);
[~, idx_min] = min(valid_results(:,6));
best_comb = valid_results(idx_min,:);

fprintf('Best Combination (minimizes time to 1000 Å):\n');
fprintf('Temperature = %.2f C \n', best_comb(1));
fprintf('Pressure = %.2f Pa\n', best_comb(2));
fprintf('H2 Concentration = %.2f SCCM\n', best_comb(3));
fprintf('Delta = %.4f m\n', best_comb(4));
fprintf('SiH4 Concentration = %.2f SCCM\n', best_comb(5));
fprintf('Time to reach 1000 Å = %.2f s\n', best_comb(6));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper Functions

function [t, thickness_A] = runSimulation(T_mod, P_mod, H2_mod, delta_mod, SiH4_mod)
    % Load baseline constants from Micron_constants1.m
    [R, ~, E_decom, A_decom, SiH4_conc, ~, tspan, ~, delta, ...
     D_SiH4, D_SiH2, D_H2, A_SiH4, E_SiH4, A_SiH2, E_SiH2, A_H2des, E_H2des, ...
     A_SiH2_growth, E_SiH2_growth, Rad, A, Thick, DS, MMs, Ratio, A_tot, ...
     A_dep_SiH4, E_dep_SiH4, A_dep_growth, E_dep_growth] = Micron_constants1();
    
    % Override baseline values with the provided parameters:
    T = T_mod;
    P = P_mod;
    H2_conc = H2_mod;
    delta = delta_mod;
    SiH4_conc = SiH4_mod;
    
    % Use an extended time span to capture slower growth scenarios
    tspan = [0, 2000];
    
    % Recalculate diffusivities (depend on T)
    D_SiH4 = (1/100^2) * (9.6e-5) * (T^1.5);
    D_SiH2 = (1/100^2) * (3.45e-5) * (T^1.5);
    D_H2   = (1/100^2) * (3.45e-5) * (T^1.5);
    
    % Initial conditions: [SiH4_g; SiH2_g; H2_g; SiH4_b; SiH2_b; H2_b]
    C0 = [SiH4_conc; 0; H2_conc; 0; 0; 0];
    
    % Reaction constant (Arrhenius relation)
    k1 = A_decom * exp(-E_decom/(R*T)) * 0.01;
    
    % Solve the ODE system
    [t, C] = ode45(@(t, C) silaneDiffusionODE(t, C, k1, R, T, P, delta, D_SiH4, D_SiH2, D_H2), tspan, C0);
    
    % Calculate deposition flux and film growth rate
k_SiH4_dep = A_dep_SiH4 * exp(-E_dep_SiH4/(R*T));
k_p_dep    = A_dep_growth * exp(-E_dep_growth/(R*T));
J = (k_p_dep * k_SiH4_dep .* C(:,4)) ./ (1 + k_SiH4_dep .* C(:,4));

% Incorporate 60% efficiency in the deposition:
growthRate_m_s = (J * 28.085 * 0.80) / 2330;
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
        % Diffusion between gas phase and boundary layer (Fick's law)
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
