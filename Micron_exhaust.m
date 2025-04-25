% SilaneGasExitCalculation.m
% This script simulates silane decomposition and reports the exit gas concentrations.

clear; clc;

%% Baseline Conditions
T     = 273 + 600;    % Temperature [K]
P     = 66.6612;      % Pressure [Pa]
H2_in = 10;           % SCCM
delta = 0.5e-3;       % m
SiH4_in = 1000;       % SCCM

%% Run Simulation
[t, thickness_A, final_gas_out] = runSilaneSimulation(T, P, H2_in, delta, SiH4_in);

%% Report Results
fprintf('Final Gas Concentrations at Reactor Outlet:\n');
fprintf('  SiH4: %.2f SCCM\n', final_gas_out(1));
fprintf('  SiH2: %.2f SCCM\n', final_gas_out(2));
fprintf('  H2  : %.2f SCCM\n', final_gas_out(3));


%% ========================================================================
function [t, thickness_A, final_gas_out] = runSilaneSimulation(T, P, H2_conc, delta, SiH4_conc)
    % Constants (simplified for this script)
    R = 8.314;  % J/mol·K
    tspan = [0, 2000];  % seconds
    
    % Kinetic constants (examples)
    A_decom = 1e6;     E_decom = 150e3;  % for SiH4 → SiH2 + H2
    A_dep_SiH4 = 1e5;  E_dep_SiH4 = 100e3;  % adsorption
    A_dep_growth = 1e4; E_dep_growth = 90e3;  % surface growth

    % Temperature-dependent diffusivities (cm^2/s → m^2/s)
    D_SiH4 = (1/100^2) * 9.6e-5 * T^1.5;
    D_SiH2 = (1/100^2) * 3.45e-5 * T^1.5;
    D_H2   = (1/100^2) * 3.45e-5 * T^1.5;

    % Initial concentrations [SiH4_g; SiH2_g; H2_g; SiH4_b; SiH2_b; H2_b]
    C0 = [SiH4_conc; 0; H2_conc; 0; 0; 0];

    % Decomposition rate constant
    k1 = A_decom * exp(-E_decom / (R * T)) * 0.01;

    % Solve ODE
    [t, C] = ode45(@(t, C) gasReactorODE(t, C, k1, R, T, P, delta, D_SiH4, D_SiH2, D_H2), tspan, C0);

    % Deposition calculation
    k_SiH4_dep = A_dep_SiH4 * exp(-E_dep_SiH4 / (R * T));
    k_p_dep    = A_dep_growth * exp(-E_dep_growth / (R * T));
    J = (k_p_dep * k_SiH4_dep .* C(:,4)) ./ (1 + k_SiH4_dep .* C(:,4));

    % Deposition rate (m/s to Å/s)
    growthRate_m_s = (J * 28.085 * 0.80) / 2330;  % g/cm³ = 2330 kg/m³
    growthRate_A_s = growthRate_m_s * 1e10;
    thickness_A = cumtrapz(t, growthRate_A_s);

    % Final gas-phase concentrations at outlet
    final_gas_out = C(end, 1:3);
end

%% ------------------------------------------------------------------------
function dCdt = gasReactorODE(~, C, k1, R, T, P, delta, D_SiH4, D_SiH2, D_H2)
    % Unpack concentrations
    SiH4_g = C(1);
    SiH2_g = C(2);
    H2_g   = C(3);
    SiH4_b = C(4);
    SiH2_b = C(5);
    H2_b   = C(6);

    % Allocate output
    dCdt = zeros(6,1);

    % Reactions in gas phase
    dCdt(1) = -k1 * SiH4_g;
    dCdt(2) =  k1 * SiH4_g;
    dCdt(3) =  2 * k1 * SiH4_g;

    % Diffusion from gas to boundary layer
    dCdt(4) = (D_SiH4 / delta) * (SiH4_g - SiH4_b);
    dCdt(5) = (D_SiH2 / delta) * (SiH2_g - SiH2_b);
    dCdt(6) = (D_H2   / delta) * (H2_g   - H2_b);
end
