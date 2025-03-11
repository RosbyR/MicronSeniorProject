function [SiH4_b, SiH2_b, H2_b] = Micron_gasphase()
% Retrieve constants
[R, T, E_decom, A_decom, SiH4_conc, H2_conc, tspan, P, delta, D_SiH4, D_SiH2, D_H2] = Micron_constants1();

% Initial conditions: [SiH4_g, Si_g, H2_g, SiH4_b, Si_b, H2_b]
C0 = [SiH4_conc; 0; H2_conc; 0; 0; 0]; 
%Decomposition Arrhenious
k1  = A_decom  * exp(-E_decom  / (R*T));
% Solve ODE system
[t, C] = ode45(@(t, C) silaneDiffusionODE(t, C, k1, R, T, P, delta, D_SiH4, D_SiH2, D_H2), tspan, C0);



% Display final concentrations in the boundary layer
disp('Final concentrations in the boundary layer:');
disp(['SiH4_b: ', num2str(C(end,4))]);
disp(['SiH2_b: ', num2str(C(end,5))]);
disp(['H2_b: ', num2str(C(end,6))]);

% Function defining the ODEs
function dCdt = silaneDiffusionODE(~, C, k1, R, T, P, delta, D_SiH4, D_SiH2, D_H2)
    % Gas Phase Concentrations
    SiH4_g = C(1);  % SiH4 in gas
    SiH2_g   = C(2);  % SiH2 in gas
    H2_g   = C(3);  % H2 in gas
    
    % Boundary Layer Concentrations
    SiH4_b = C(4);  % SiH4 in boundary layer
    SiH2_b   = C(5);  % SiH2 in boundary layer
    H2_b   = C(6);  % H2 in boundary layer

    % Define ODEs
    dCdt = zeros(6,1);

    % Gas Phase Reactions
    dCdt(1) = -k1 * SiH4_g;      % d[SiH4_g]/dt = -k1 * [SiH4_g]
    dCdt(2) = k1 * SiH4_g;       % d[SiH2_g]/dt = +k1 * [SiH4_g]
    dCdt(3) = 2 * k1 * SiH4_g;   % d[H2_g]/dt = +2 * k1 * [SiH4_g]

    % Diffusion into Boundary Layer (Fick's Law)   
  
    dCdt(4) = (D_SiH4 / delta) * (SiH4_g - SiH4_b) * (P / (R * T));  % SiH4_b
    dCdt(5) = (D_SiH2 / delta) * (SiH2_g - SiH2_b) * (P / (R * T));        % SiH2_b
    dCdt(6) = (D_H2 / delta) * (H2_g - H2_b) * (P / (R * T));        % H2_b

end
end