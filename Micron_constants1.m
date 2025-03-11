function[R, T, E_decom, A_decom, SiH4_conc, H2_conc, tspan, P, delta, D_SiH4, D_SiH2, D_H2, ...
 A_SiH4, E_SiH4, A_SiH2, E_SiH2, A_H2des, E_H2des, A_SiH2_growth, E_SiH2_growth, Rad, A, Thick,...
 DS, MMs, Ratio] = Micron_constants1();

%Assumptions
%Uniform Heat Distribution, Ideal Gas Law, Surface Limited Diffusion,
%Uniform Film Thickness, Constant Boundary Layer Thickness, Constant
%Temperature, Constant Pressure, Gas Phase Concentration >> Boundary Layer,
%Assuming Fick's First Law for diffusion, 





% Universal Gas Constant
    R = 8.314;       % J/(mol*K) **
    
    % System Conditions
    T = 273 + 600;        % Temperature (K)
    P = 66.6612;      % Pressure (Pa)
    delta = 1e-1;    % Boundary layer thickness (m)  
    
    % Gas Phase Constants    Journal of Chemical Physics
    %Gordon, M. S., Gano, D. R., Binkley, J. S., & Frisch, M. J. (1985). Thermal decomposition of silane. Journal of the American Chemical Society, 107(23)
    E_decom = 209000; %Activiation Energy
    A_decom = 1.2e15; %Pre-exponential Factor
    
    
    % Initial Concentrations 
    SiH4_conc = .10; %mol/m^3
    H2_conc = 10; %mol/m^3
    
    % Time Span
    tspan = [0, 200]; % Time span for reaction (seconds)

    % Diffusivities 
    D_SiH4 = 4; %boundary layer
    D_SiH2 = 5;
    D_H2   = 1;

    % Surface adsorption constants: DOI:10.1088/0268-1242/6/4/009 
    A_SiH4 = 38; %K in units of kmol^-1 m^3
    E_SiH4 = 58e6; 
    A_SiH2 = 38;
    E_SiH2 = 58e6;
    A_H2des = 38;
    E_H2des = 47e6;




    %Growth Constants DOI:10.1088/0268-1242/6/4/009 
    %K in Units of kmol m^-2 s^-1
    A_SiH2_growth = 69; %
    E_SiH2_growth = -158; %

    Rad = 150; %Wafer Radius in MM
    A = pi * Rad^2; %Wafer Area in MM^2
    Thick = 1000; %Film Thickness in Angstroms
    DS = 2330; % Silane Density in Kg/M^3
    MMs = 28.085; % Molar Mass Silane in g/Mol
    Ratio = H2_conc / SiH4_conc;


end
