%% Micron_deposition_best_case.m
% Runs the deposition simulation using the best-case parameters
% and generates all relevant plots and output.

function [t, thickness_A, exhaust] = runSimulation()
    % Load constants and best-case parameters
    [R, T, E_decom, A_decom, SiH4_conc, H2_conc, ~, P, delta, D_SiH4, D_SiH2, D_H2, ...
     A_SiH4, E_SiH4, A_SiH2, E_SiH2, A_H2des, E_H2des, A_SiH2_growth, E_SiH2_growth, Rad, A, Thick,...
     DS, MMs, Ratio, A_tot, A_dep_SiH4, E_dep_SiH4, A_dep_growth, E_dep_growth, ...
     SiH4_sccm, V_chamber, F_in] = Micron_constants1();

    % Override for best-case scenario
    T = 725;              % K
    P = 1;                % Pa
    delta = 0.005;        % m
    H2_conc = 100;        % mol/m^3
    SiH4_conc = 1000;     % mol/m^3

    tspan = [0, 500];    % Extended time span for slow growth

    % Recalculate diffusivities
    D_SiH4 = (1/100^2) * (9.6e-5) * (T^1.5);
    D_SiH2 = (1/100^2) * (3.45e-5) * (T^1.5);
    D_H2   = (1/100^2) * (3.45e-5) * (T^1.5);

    % Initial conditions: [SiH4_g; SiH2_g; H2_g; SiH4_b; SiH2_b; H2_b]
    C0 = [SiH4_conc; 0; H2_conc; 0; 0; 0];

    % Reaction constant
    k1 = A_decom * exp(-E_decom/(R*T)) * 0.01;

    % Solve the ODE system
    [t, C] = ode45(@(t, C) silaneDiffusionODE(t, C, k1, delta, D_SiH4, D_SiH2, D_H2), tspan, C0);

    % Deposition kinetics
    k_SiH4_dep = A_dep_SiH4 * exp(-E_dep_SiH4/(R*T));
    k_p_dep    = A_dep_growth * exp(-E_dep_growth/(R*T));
    J = (k_p_dep * k_SiH4_dep .* C(:,4)) ./ (1 + k_SiH4_dep .* C(:,4));

    % Growth rate and film thickness
    growthRate_m_s = (J * 28.085 * 0.80) / 2330;
    growthRate_A_s = growthRate_m_s * 1e10;
    thickness_A = cumtrapz(t, growthRate_A_s);

    %% Plot Results
    figure;
    plot(t, thickness_A, 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Film Thickness (Angstrom)');
    title('Film Growth Over Time'); grid on;

    

    figure;
    plot(t, C(:,4), t, C(:,5), t, C(:,6), 'LineWidth', 2);
    legend('SiH_4', 'SiH_2', 'H_2');
    title('Boundary Layer Species Concentrations'); xlabel('Time (s)'); ylabel('mol/m^3'); grid on;

    figure;
    plot(t, J, 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Deposition Flux (kmol/m^2/s)');
    title('Deposition Flux Over Time'); grid on;

    %% Exhaust Composition
    C_out = [C(end,1); C(end,2); C(end,3)];
    total_C = sum(C_out);
    mole_frac = C_out / total_C;
    F_total = 1e-3;
    molar_outflow = mole_frac * F_total;

    exhaust.H2   = molar_outflow(3);
    exhaust.SiH4 = molar_outflow(1);
    exhaust.SiH2 = molar_outflow(2);

    exhaustTable = table({'H2'; 'SiH4'; 'SiH2'}, ...
                         [mole_frac(3); mole_frac(1); mole_frac(2)], ...
                         [exhaust.H2; exhaust.SiH4; exhaust.SiH2], ...
                         'VariableNames', {'Species', 'MoleFraction', 'MolarFlowrate'});
    disp('--- Exhaust Composition ---');
    disp(exhaustTable);

    %% Nested ODE
    function dCdt = silaneDiffusionODE(~, C, k1, delta, D1, D2, D3)
        dCdt = zeros(6,1);
        dCdt(1) = -k1 * C(1);
        dCdt(2) =  k1 * C(1);
        dCdt(3) = 2 * k1 * C(1);
        dCdt(4) = (D1 / delta) * (C(1) - C(4));
        dCdt(5) = (D2 / delta) * (C(2) - C(5));
        dCdt(6) = (D3 / delta) * (C(3) - C(6));
    end
end