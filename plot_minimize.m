% Temperature Sensitivity Analysis Script
clear; clc;

T_range = 700:1:750;  % Temperature range to sweep
thicknesses = zeros(size(T_range));

% Get reference thickness at T_ref = 725 K
T_ref = 725;
[~, thickness_ref, ~] = runSimulationAtTemp(T_ref);
final_ref = thickness_ref(end);

% Sweep temperatures and compute final thickness
for i = 1:length(T_range)
    T = T_range(i);
    [~, thickness, ~] = runSimulationAtTemp(T);
    thicknesses(i) = thickness(end);
end

% Calculate percent variation from reference
percent_diff = 100 * (thicknesses - final_ref) / final_ref;

% Find indices where variation is within ±3%
within_range_idx = find(abs(percent_diff) <= 1.5);
valid_temperatures = T_range(within_range_idx);

fprintf('Temperature range with ±3%% thickness variation: %d K to %d K\n', ...
        min(valid_temperatures), max(valid_temperatures));

% Plot
figure;
plot(T_range, percent_diff, '-o', 'LineWidth', 2);
yline(3, '--r'); yline(-3, '--r');
xlabel('Temperature (C)');
ylabel('Percent Thickness Change');
title('Temperature Sensitivity of Film Thickness');
grid on;

% --- Local Function ---
function [t, thickness_A, exhaust] = runSimulationAtTemp(T)
    % Load constants and override temperature
    [R, ~, E_decom, A_decom, SiH4_conc, H2_conc, ~, P, delta, D_SiH4, D_SiH2, D_H2, ...
     A_SiH4, E_SiH4, A_SiH2, E_SiH2, A_H2des, E_H2des, A_SiH2_growth, E_SiH2_growth, Rad, A, Thick,...
     DS, MMs, Ratio, A_tot, A_dep_SiH4, E_dep_SiH4, A_dep_growth, E_dep_growth, ...
     SiH4_sccm, V_chamber, F_in] = Micron_constants1();

    delta = 0.005; P = 1;
    H2_conc = 100; SiH4_conc = 1000;

    tspan = [0, 500];
    D_SiH4 = (1/100^2) * (9.6e-5) * (T^1.5);
    D_SiH2 = (1/100^2) * (3.45e-5) * (T^1.5);
    D_H2   = (1/100^2) * (3.45e-5) * (T^1.5);

    C0 = [SiH4_conc; 0; H2_conc; 0; 0; 0];
    k1 = A_decom * exp(-E_decom/(R*T)) * 0.01;

    [t, C] = ode45(@(t, C) ode(t, C, k1, delta, D_SiH4, D_SiH2, D_H2), tspan, C0);

    k_SiH4_dep = A_dep_SiH4 * exp(-E_dep_SiH4/(R*T));
    k_p_dep    = A_dep_growth * exp(-E_dep_growth/(R*T));
    J = (k_p_dep * k_SiH4_dep .* C(:,4)) ./ (1 + k_SiH4_dep .* C(:,4));

    growthRate_m_s = (J * 28.085 * 0.80) / 2330;
    growthRate_A_s = growthRate_m_s * 1e10;
    thickness_A = cumtrapz(t, growthRate_A_s);

    exhaust = [];  % Not used here

    function dCdt = ode(~, C, k1, delta, D1, D2, D3)
        dCdt = zeros(6,1);
        dCdt(1) = -k1 * C(1);
        dCdt(2) =  k1 * C(1);
        dCdt(3) = 2 * k1 * C(1);
        dCdt(4) = (D1 / delta) * (C(1) - C(4));
        dCdt(5) = (D2 / delta) * (C(2) - C(5));
        dCdt(6) = (D3 / delta) * (C(3) - C(6));
    end
end
