function [R, T, k1, SiH4_conc, tspan, P, delta, D_SiH4, D_SiH2, ...
    D_H2,A_SiH4, E_SiH4, A_SiH2, E_SiH2, A_H2des, E_H2des] = Micron_constants1()
    % Universal Gas Constant
    R = 8.314;       % J/(mol*K)
    
    % System Conditions
    T = 273 + 600;        % Temperature (K)
    P = 66.6612;      % Pressure (Pa)
    delta = 1e-4;    % Boundary layer thickness (m)  
    %**No Idea if this is a reasonable boundary layer
    
    % Reaction Rate Constants
    k1 = 1e5;        % Rate constant for SiH4 decomposition (1/s)

    % Initial Concentrations
    SiH4_conc = 1;   % Initial concentration of SiH4 (mol/m^3)
    
    % Time Span
    tspan = [0 1e-4]; % Time span for reaction (seconds)

    % **Dummy Diffusivities** (m²/s) - Example values for gases
    D_SiH4 = 1e-5;  % Diffusivity of SiH4 (assumed)
    D_SiH2   = 5e-6;  % Diffusivity of SiH2 (assumed)
    D_H2   = 7e-5;  % Diffusivity of H2 (assumed)


    %Surface adsorption constants
    A_SiH4 = 8.39e26;
    E_SiH4 = 37450.0;
    A_SiH2 = 8.39e27;
    E_SiH2 = 37450.0;
    A_H2des = 1.75e20;
    E_H2des = 47000.0;

end