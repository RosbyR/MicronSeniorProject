function[R, T, k1, SiH4_conc, tspan, P, delta, D_SiH4, D_SiH2, D_H2, ...
 A_SiH4, E_SiH4, A_SiH2, E_SiH2, A_H2des, E_H2des] = Micron_constants1();

    
% Universal Gas Constant
    R = 8.314;       % J/(mol*K)
    
    % System Conditions
    T = 273 + 600;        % Temperature (K)
    P = 66.6612;      % Pressure (Pa)
    delta = 1e-4;    % Boundary layer thickness (m)  
    
    % Reaction Rate Constants (unused in your main code)
    k1 = 1e5;        
    
    % Initial Concentrations (unused in your main code)
    SiH4_conc = 1;   
    
    % Time Span
    tspan = [0, 1e-4]; % Time span for reaction (seconds)

    % Diffusivities (unused in your main code)
    D_SiH4 = 1e-5;
    D_SiH2 = 5e-6;
    D_H2   = 7e-5;

    % Surface adsorption constants
    A_SiH4 = 8.39e26;
    E_SiH4 = 37450.0;
    A_SiH2 = 8.39e27;
    E_SiH2 = 37450.0;
    A_H2des = 1.75e20;
    E_H2des = 47000.0;
end
