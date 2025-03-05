function[R, T, k1, SiH4_conc, tspan, P, delta, D_SiH4, D_SiH2, D_H2, ...
 A_SiH4, E_SiH4, A_SiH2, E_SiH2, A_H2des, E_H2des] = Micron_constants1();

%Assumptions
%Uniform Heat Distribution,



% Universal Gas Constant
    R = 8.314;       % J/(mol*K) **
    
    % System Conditions
    T = 273 + 600;        % Temperature (K)
    P = 66.6612;      % Pressure (Pa)
    delta = 1e-3;    % Boundary layer thickness (m)  
    
    % Gas Phase Constants    Journal of Chemical Physics
    %Gordon, M. S., Gano, D. R., Binkley, J. S., & Frisch, M. J. (1985). Thermal decomposition of silane. Journal of the American Chemical Society, 107(23)
    E_decom = 209000; %Activiation Energy
    A_decom = 1.2e15; %Pre-exponential Factor
    k1 = 1e5;        %Silane Decomposition
    
    % Initial Concentrations 
    SiH4_conc = .1; %mol/m^3
    
    % Time Span
    tspan = [0, 30]; % Time span for reaction (seconds)

    % Diffusivities 
    D_SiH4 = 1e-2; %boundary layer
    D_SiH2 = 5e-2;
    D_H2   = 7e-2;

    % Surface adsorption constants
    A_SiH4 = 8.39e26; %paper
    E_SiH4 = 37450.0; 
    A_SiH2 = 8.39e27;
    E_SiH2 = 37450.0;
    A_H2des = 1.75e20;
    E_H2des = 47000.0;
end
