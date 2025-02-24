function [R, T, k1, SiH4_conc, tspan, P, delta, D_SiH4, D_Si, D_H2] = Micron_constants()
    % Universal Gas Constant
    R = 8.314;       % J/(mol*K)
    
    % System Conditions
    T = 273 + 600;        % Temperature (K)
    P = 66.6612;      % Pressure (Pa)
    delta = 1e-4;    % Boundary layer thickness (m)
    
    % Reaction Rate Constants
    k1 = 1e5;        % Rate constant for SiH4 decomposition (1/s)

    % Initial Concentrations
    SiH4_conc = 1;   % Initial concentration of SiH4 (mol/m^3)
    
    % Time Span
    tspan = [0 1e-4]; % Time span for reaction (seconds)

    % **Dummy Diffusivities** (m²/s) - Example values for gases
    D_SiH4 = 1e-5;  % Diffusivity of SiH4 (assumed)
    D_Si   = 5e-6;  % Diffusivity of Si (assumed)
    D_H2   = 7e-5;  % Diffusivity of H2 (assumed)
end