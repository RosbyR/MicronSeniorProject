%% Revised Model: Solid Deposition with H2 Dynamics, Diffusion, and H2 Inhibition

% Define constants and parameters
R       = 8.314;                   % Universal gas constant (J/(mol*K))
T       = 600 + 273.15;            % Temperature (K)
gamma   = 0.5;                     % Surface tension (J/m^2)
V_bar   = 2.5e-29;                 % Atomic volume (m^3)

% Kinetic constants
k1      = 1e5;                     % Rate constant for SiH4 -> SiH2 (1/s)
k2      = 1e3;                     % Rate constant for SiH2 consumption (1/s)
k_ads   = 1e1;                     % Nucleation (initiation) rate constant
k_des   = 1e-5;                    % Desorption rate constant for solid silicon (1/s)
k_growth= 1e-3;                    % Growth rate constant on existing nuclei

% Assume excess silane (constant concentration, mol/m^3)
SiH4_conc = 1e-2;  
% Steady-state SiH2 concentration (determined by the reaction constants)
SiH2_const = (k1/k2) * SiH4_conc;

% Diffusion parameters for H2
D0      = 1e-5;                   % Pre-exponential factor for diffusivity (m^2/s)
Ea_D    = 45e3;                   % Activation energy for diffusivity (J/mol)
delta   = 1e-5;                   % Boundary layer thickness (m)
DH2     = D0 * exp(-Ea_D/(R*T));  % Diffusivity for H2 (m^2/s)
k_m_H2  = DH2 / delta;            % Mass transfer coefficient for H2 (m/s)

% Inhibition parameter for H2 (adjust to control how strongly H2 blocks growth)
alpha_inh = 1e-3;  % Units: m^3/mol

% Initial conditions:
% y(1) = [Si_s] (solid silicon) - start with a tiny nucleation seed.
% y(2) = [H2] - assume initially zero.
Si_s0 = 1e-6;    
H2_0  = 0;       
y0 = [Si_s0; H2_0];

% Time span (seconds)
tspan = [0, 1000];

% ODE solver options
options = odeset('NonNegative', 1, 'RelTol',1e-6, 'AbsTol',1e-9);

% Solve the ODE system
[t, y] = ode15s(@(t,y) cvd_kinetics_H2(t, y, k1, k2, k_ads, k_des, k_growth, ...
    SiH4_conc, SiH2_const, k_m_H2, gamma, V_bar, R, T, alpha_inh), ...
    tspan, y0, options);

% Plot the results
figure;
plot(t, y(:,1), 'k', 'LineWidth', 2, 'DisplayName', '[Si_s] (Solid Silicon)');
hold on;
%plot(t, y(:,2), 'b', 'LineWidth', 2, 'DisplayName', '[H_2]');
%plot(t, SiH4_conc*ones(size(t)), 'r--', 'LineWidth', 1.5, 'DisplayName', '[SiH_4] (Excess)');
xlabel('Time (s)');
ylabel('Concentration (mol/m^3)');
legend;
title('Solid Silicon Deposition with H_2 Dynamics and H_2 Inhibition');
hold off;

function dydt = cvd_kinetics_H2(t, y, k1, k2, k_ads, k_des, k_growth, ...
                                SiH4_conc, SiH2_const, k_m_H2, gamma, V_bar, R, T, alpha_inh)
    % y(1) = [Si_s] (solid silicon)
    % y(2) = [H2]
    
    % --- Nucleation Rate Calculation ---
    % Assume a reference equilibrium concentration/partial pressure p_eq:
    p_eq = 1e-3;  
    supersat = SiH4_conc / p_eq;
    if (supersat <= 1) || isnan(log(supersat))
        nucleation_rate = 0;
    else
        r_star = (2 * gamma * V_bar) / (R * T * log(supersat));
        nucleation_rate = k_ads * r_star;
    end
    
    % --- Inhibition Factor due to H2 ---
    % High H2 concentration reduces the effective growth rate.
    inhibition_factor = 1 / (1 + alpha_inh * y(2));
    
    % --- Solid Silicon (Si_s) Deposition ---
    % Growth is driven by the available SiH2 and existing Si_s but is inhibited by H2.
    dSi_s_dt = nucleation_rate + k_growth * SiH2_const * y(1) * inhibition_factor;
    
    % --- H2 Dynamics ---
    % Assume that for each conversion of SiH4 to SiH2, H2 is produced.
    % Production term: k1 * [SiH4] (constant since SiH4 is in excess).
    % Loss: Diffusion through the boundary layer, modeled as k_m_H2 * [H2].
    dH2_dt = k1 * SiH4_conc - k_m_H2 * y(2);
    
    dydt = [dSi_s_dt; dH2_dt];
end

