
function[R, T, E_decom, A_decom, SiH4_conc, H2_conc, tspan, P, delta, D_SiH4, D_SiH2, D_H2, ...
 A_SiH4, E_SiH4, A_SiH2, E_SiH2, A_H2des, E_H2des, A_SiH2_growth, E_SiH2_growth, Rad, A, Thick,...
 DS, MMs, Ratio, A_tot, A_dep_SiH4, E_dep_SiH4, A_dep_growth, E_dep_growth, ...
 SiH4_sccm, V_chamber, F_in] = Micron_constants1()



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Assumptions
%Ideal Gas, Uniform boundary layer, uniform temperature, constant pressure
%Steady State except for surface growth, 20% Loss to walls of reactor,
%Surface Limited Diffusion, 
    % Universal Gas Constant
    R = 8.314;       % J/(mol*K)
    
    % System Conditions
    T = 273 + 600;        % Temperature (K)
    P = 1;          % Pressure (Pa)
    delta = .005;       % Boundary layer thickness (m)
    
    % Gas Phase Constants
    E_decom = 209000;     % Activation Energy (J/mol)
    A_decom = 1.2e15;     % Pre-exponential Factor (1/s)
    
    % Initial Concentrations 
    SiH4_conc = 1000;   % mol/m^3
    H2_conc = 100;        % mol/m^3
    
    % Time Span
    tspan = [0, 600];     % Time span for reaction (s)
    
    % Diffusivities 
    D_SiH4 = (1/100^2)*(9.6e-5)*(T^1.5); % m^2/s
    D_SiH2 = (1/100^2)*(3.45e-5)*(T^1.5);
    D_H2   = (1/100^2)*(3.45e-5)*(T^1.5);
    
    % Surface adsorption constants
    A_SiH4 = 38; E_SiH4 = 58e6;
    A_SiH2 = 38; E_SiH2 = 58e6;
    A_H2des = 38; E_H2des = 47e6;
    
    % Growth Constants
    A_SiH2_growth = 69;
    E_SiH2_growth = -158;
    
    Rad = 150; % mm
    A = pi * Rad^2; % mm^2
    Thick = 1000; % Å
    DS = 2330; % kg/m^3
    MMs = 28.085; % g/mol
    Ratio = H2_conc / SiH4_conc;
    A_m2 = A / (1000 * 1000); % m^2
    A_tot = A_m2 * 2 * 100; % total wafer area
    
    A_dep_SiH4   = 50;
    E_dep_SiH4   = 45e3;
    A_dep_growth = 80;
    E_dep_growth = 130e3;

    % Inflow Settings
    SiH4_sccm = 1000;                % Flow in sccm
    V_chamber = 0.1;                 % Chamber volume in m^3
    mol_per_sccm = 7.45e-7;          % mol/s per sccm at STP
    F_in = SiH4_sccm * mol_per_sccm; % mol/s

end
