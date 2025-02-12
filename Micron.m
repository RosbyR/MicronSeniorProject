% Define constants and new parameters
A = 1.5e10;      % Pre-exponential factor in m/s
Ea = 53.7e3;     % Activation energy in J/mol
R = 8.314;       % Universal gas constant in J/mol·K
T = 600 + 273.15; % Temperature in Kelvin
P = 20;           % Increased Pressure to increase supersaturation
gamma = 0.5;      % Surface tension (example value, in J/m^2)
V_bar = 2.5e-29;  % Atomic volume (example value, in m^3)

% Rate constants for reactions
k1 = 1e5;   % Rate constant for SiH4 -> SiH2
k2 = 1e3;   % Rate constant for SiH2 decomposition
k_ads = 1e1; % Increased adsorption rate constant (much faster deposition)
k_des = 1e-5; % Reduced desorption rate constant (slower loss)
k_H2 = 1e-3;  % Hydrogen desorption rate constant

% Temperature-dependent diffusivity
D0 = 1e-5;   % Pre-exponential factor for diffusivity (m^2/s)
Ea_D = 45e3;  % Activation energy for diffusivity (J/mol)

% Diffusivity of components (m^2/s)
DS4 = D0 * exp(-Ea_D / (R * T));  
DS2 = D0 * exp(-Ea_D / (R * T));  
DH2 = D0 * exp(-Ea_D / (R * T));  

% Boundary layer thickness (meters)
delta = 1e-5;

% Mass transfer coefficients (m/s)
k_m_SiH4 = DS4 / delta;  % Mass transfer coefficient for SiH4
k_m_SiH2 = DS2 / delta;  % Mass transfer coefficient for SiH2

% Define the system of differential equations for CVD kinetics
function dydt = cvd_kinetics(t, y, A, Ea, R, T, k1, k2, k_ads, k_des, k_H2, k_m_SiH4, k_m_SiH2, gamma, V_bar, SiH4_influx)
    % Extract concentrations
    SiH4 = max(y(1), 1e-12);  % Silane concentration, preventing zero
    SiH2 = max(y(2), 1e-12);  % Silylene concentration
    H2 = max(y(3), 1e-12);    % Hydrogen concentration
    Si_s = max(y(4), 1e-6);      % Start with a larger initial surface silicon concentration

    % Arrhenius rate constant
    k = A * exp(-Ea / (R * T));

    % Mass transfer limited flux (surface transport)
    J_SiH4 = k_m_SiH4 * (SiH4 - Si_s); % Flux of SiH4
    J_SiH2 = k_m_SiH2 * (SiH2 - Si_s); % Flux of SiH2

    % Nucleation rate (based on critical cluster size equation)
    p_eq = 1e-3;  % Equilibrium pressure (example value in atm)
    supersaturation = SiH4 / p_eq;  % Supersaturation (example)
    r_star = (2 * gamma * V_bar) / (R * T * log(supersaturation));  % Critical cluster size
    nucleation_rate = k_ads * r_star;  % Nucleation rate is proportional to the cluster size

    % Differential equations
    dydt = zeros(4, 1);
    
    % Constant influx of SiH4 into the system (no depletion)
    dydt(1) = SiH4_influx - J_SiH4 - k * SiH4;  % SiH4 decay due to reaction + transport, and constant influx
    
    dydt(2) = k1 * SiH4 - k2 * SiH2 - J_SiH2;    % SiH2 formation and decomposition
    dydt(3) = k_H2 * H2 - k_H2 * H2;              % H2 production and desorption (excess hydrogen desorbs)
    dydt(4) = nucleation_rate * Si_s - k_des * Si_s;     % Surface silicon deposition with nucleation term
end

% Initial conditions - Adjusted to realistic values
y0 = [1.0; 1e-3; 1e-6; 0.10];  % Start SiH4 at 1 mol/m^3, SiH2 at 1e-3 mol/m^3, Si_s at 1 mol/m^3 (initial solid silicon)

% Constant SiH4 influx (mol/m^3/s)
SiH4_influx = 1e-3;  % Influx rate of SiH4 into the system (mol/m^3/s)

% Solve the ODE using a stiff solver
options = odeset('NonNegative', [1 2 3 4], 'RelTol', 1e-6, 'AbsTol', 1e-9);
[t, y] = ode15s(@(t, y) cvd_kinetics(t, y, A, Ea, R, T, k1, k2, k_ads, k_des, k_H2, k_m_SiH4, k_m_SiH2, gamma, V_bar, SiH4_influx), [0, 1000], y0, options);

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
