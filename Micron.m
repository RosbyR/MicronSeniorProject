function dydt = cvd_kinetics(t, y, A, Ea, R, T, k_ads, k_des, k_H2)
    SiH4 = y(1);  % Silane concentration
    SiH2 = y(2);  % Silylene concentration
    H2 = y(3);    % Hydrogen concentration
    Si_s = y(4);  % Surface silicon concentration

    % Rate constant from Arrhenius equation
    k = A * exp(-Ea / (R * T));  

    % Differential equations for each species
    dydt = zeros(4, 1);
    dydt(1) = -k * SiH4;                      % Silane concentration decay
    dydt(2) = k1 * SiH4 - k2 * SiH2;           % Silylene dynamics
    dydt(3) = k_H2 * H2;                      % Hydrogen desorption
    dydt(4) = k_ads * SiH4 - k_des * Si_s;    % Surface silicon dynamics
end

% Initial concentrations and parameters
y0 = [1.0; 0; 0; 0];  % Initial concentrations [SiH4, SiH2, H2, Si_s]
[t, y] = ode45(@(t, y) cvd_kinetics(t, y, A, Ea, R, T, k_ads, k_des, k_H2), [0, 1000], y0);

% Plot results
figure;
plot(t, y(:, 1), 'r', 'DisplayName', '[SiH4]');
hold on;
plot(t, y(:, 2), 'g', 'DisplayName', '[SiH2]');
plot(t, y(:, 3), 'b', 'DisplayName', '[H2]');
plot(t, y(:, 4), 'k', 'DisplayName', '[Si_s]');
xlabel('Time (s)');
ylabel('Concentration (mol/m^3)');
legend;
