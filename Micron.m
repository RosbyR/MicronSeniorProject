% Define constants
A = 1.5e10;      % Pre-exponential factor in m/s
Ea = 53.7e3;     % Activation energy in J/mol
R = 8.314;       % Universal gas constant in J/mol·K
T = 600 + 273.15; % Temperature in Kelvin

% Rate constants
k1 = 1e5;   % Rate constant for SiH4 -> SiH2
k2 = 1e3;   % Rate constant for SiH2 decomposition
k_ads = 1e-3; % Adsorption rate constant
k_des = 1e-3; % Desorption rate constant
k_H2 = 1e-3;  % Hydrogen desorption rate constant

% Diffusivity of components (m^2/s)
DS4 = 1e-5;  
DS2 = 1e-5;  
DH2 = 1e-5;  

% Boundary layer thickness (meters)
delta = 1e-5;

% Mass transfer coefficients (m/s)
k_m_SiH4 = DS4 / delta;  % Mass transfer coefficient for SiH4
k_m_SiH2 = DS2 / delta;  % Mass transfer coefficient for SiH2

% Define the system of differential equations for CVD kinetics
function dydt = cvd_kinetics(t, y, A, Ea, R, T, k1, k2, k_ads, k_des, k_H2, k_m_SiH4, k_m_SiH2)
    % Extract concentrations
    SiH4 = max(y(1), 1e-12);  % Silane concentration, preventing zero
    SiH2 = max(y(2), 1e-12);  % Silylene concentration
    H2 = max(y(3), 1e-12);    % Hydrogen concentration
    Si_s = max(y(4), 0);      % Surface silicon concentration

    % Arrhenius rate constant
    k = A * exp(-Ea / (R * T));

    % Mass transfer limited flux (surface transport)
    J_SiH4 = k_m_SiH4 * (SiH4 - Si_s); % Flux of SiH4
    J_SiH2 = k_m_SiH2 * (SiH2 - Si_s); % Flux of SiH2

    % Differential equations
    dydt = zeros(4, 1);
    dydt(1) = -J_SiH4 - k * SiH4;                % SiH4 decay due to reaction + transport
    dydt(2) = k1 * SiH4 - k2 * SiH2 - J_SiH2;    % SiH2 formation and decomposition
    dydt(3) = k_H2 * H2;                         % H2 production (no loss term)
    dydt(4) = k_ads * J_SiH4 - k_des * Si_s;     % Surface silicon deposition
end

% Initial conditions
y0 = [1.0; 1e-6; 0; 0];  % Start SiH2 at a small value for numerical stability

% Solve the ODE using a stiff solver
options = odeset('NonNegative', [1 2 3 4], 'RelTol', 1e-6, 'AbsTol', 1e-9);
[t, y] = ode15s(@(t, y) cvd_kinetics(t, y, A, Ea, R, T, k1, k2, k_ads, k_des, k_H2, k_m_SiH4, k_m_SiH2), [0, 1000], y0, options);

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
hold off;
