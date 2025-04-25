% DepositionSensitivitySweep.m
% Sensitivity analysis on deposition rate vs. temperature, SiH4 flow, and H2:SiH4 ratio

clear; clc;

% Define input ranges
temperatures = [673, 773, 1073];           % in Kelvin
sih4_flows   = [500, 1000, 1500];          % SCCM
h2_ratios    = [0.0, 0.5, 1.0, 0.1];       % H2:SiH4 ratio

% Allocate results
results = [];

% Loop over each combination
for iT = 1:length(temperatures)
    for iF = 1:length(sih4_flows)
        for iR = 1:length(h2_ratios)
            
            % Set inputs
            T = temperatures(iT);
            SiH4_sccm = sih4_flows(iF);
            H2_sccm   = h2_ratios(iR) * SiH4_sccm;
            delta     = 0.005;         % Default value for boundary layer thickness (m)
            P         = 1;         % Pressure in Pa
            
            try
                % Call your deposition model
                [t, thickness_A] = runSimulation(T, P, H2_sccm, delta, SiH4_sccm);
                
                % Compute deposition rate as average slope (Å/s)
                if length(t) >= 2
                    dep_rate = (thickness_A(end) - thickness_A(1)) / (t(end) - t(1));
                else
                    dep_rate = NaN;
                end
            catch
                dep_rate = NaN;
            end

            % Store result
            results = [results; T, SiH4_sccm, H2_sccm, h2_ratios(iR), dep_rate];
        end
    end
end

% Create table
dep_table = array2table(results, ...
    'VariableNames', {'Temperature_K', 'SiH4_SCCM', 'H2_SCCM', 'H2_SiH4_Ratio', 'DepRate_A_per_s'});

disp(dep_table);

% Optional: Plotting
figure;
gscatter(dep_table.H2_SiH4_Ratio, dep_table.DepRate_A_per_s, dep_table.Temperature_K);
xlabel('H2:SiH4 Ratio');
ylabel('Deposition Rate (Å/s)');
title('Deposition Rate Sensitivity to H2:SiH4 Ratio');
grid on;
legend('Location','best');


function [t, thickness_A] = runSimulation(T, P, H2_conc, delta, SiH4_conc)

    % Load constants
    [R, ~, E_decom, A_decom, ~, ~, ~, ~, ~, ...
     D_SiH4, D_SiH2, D_H2, A_SiH4, E_SiH4, A_SiH2, E_SiH2, ...
     A_H2des, E_H2des, A_SiH2_growth, E_SiH2_growth, Rad, A, Thick, ...
     DS, MMs, Ratio, A_tot, A_dep_SiH4, E_dep_SiH4, ...
     A_dep_growth, E_dep_growth] = Micron_constants1();

    % Time span
    tspan = [0 2000];

    % Recalculate diffusivities
    D_SiH4 = (1/100^2) * 9.6e-5 * T^1.5;
    D_SiH2 = (1/100^2) * 3.45e-5 * T^1.5;
    D_H2   = (1/100^2) * 3.45e-5 * T^1.5;

    % Initial concentrations: [SiH4_g, SiH2_g, H2_g, SiH4_b, SiH2_b, H2_b]
    C0 = [SiH4_conc; 0; H2_conc; 0; 0; 0];

    % Arrhenius rate constant
    k1 = A_decom * exp(-E_decom / (R*T)) * 0.01;

    % Solve ODE
    [t, C] = ode45(@(t, C) silaneODE(t, C, k1, delta, D_SiH4, D_SiH2, D_H2), tspan, C0);

    % Calculate deposition flux
    k_SiH4_dep = A_dep_SiH4 * exp(-E_dep_SiH4 / (R*T));
    k_p_dep    = A_dep_growth * exp(-E_dep_growth / (R*T));
    J = (k_p_dep * k_SiH4_dep .* C(:,4)) ./ (1 + k_SiH4_dep .* C(:,4));

    % Growth rate with 80% efficiency
    growthRate_m_s = (J * 28.085 * 0.80) / 2330;
    growthRate_A_s = growthRate_m_s * 1e10;  % to Å/s

    % Integrate for thickness
    thickness_A = cumtrapz(t, growthRate_A_s);

    % --- Nested ODE Function ---
    function dCdt = silaneODE(~, C, k1, delta, D1, D2, D3)
        dCdt = zeros(6,1);
        dCdt(1) = -k1 * C(1);               % SiH4_g
        dCdt(2) =  k1 * C(1);               % SiH2_g
        dCdt(3) =  2 * k1 * C(1);           % H2_g
        dCdt(4) = (D1/delta) * (C(1)-C(4)); % SiH4_b
        dCdt(5) = (D2/delta) * (C(2)-C(5)); % SiH2_b
        dCdt(6) = (D3/delta) * (C(3)-C(6)); % H2_b
    end
end
