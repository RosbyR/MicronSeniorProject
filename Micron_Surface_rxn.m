% Retrieve constants
[R, T, A_SiH4, E_SiH4, A_SiH2, E_SiH2, A_H2des, E_H2des] = Micron_constants1();


[SiH4_b, SiH2_b, H2_b] = Micron_gasphase();

%Compute rate based on first order kinetics
k_silane   = A_SiH4   * exp(-E_SiH4   / (R*T));
k_disilane = A_SiH2 * exp(-E_SiH2 / (R*T));
k_H2des    = A_H2des    * exp(-E_H2des    / (R*T));



y0 = [SiH4_b; SiH2_b; H2_b; 0;]

function dydt = SurfaceRxnODE(~, y, R, T, k_silane, k_disilane, k_H2des)
    % Unpack the variables (example: partial pressures or concentrations)
    p_SiH4  = y(1);
    p_Si2H6 = y(2);
    p_H2    = y(3);
    theta_H = y(4);  % surface coverage, for instance

% Example ODEs (placeholders):
    dp_SiH4_dt  = -k_silane * p_SiH4;  
    dp_Si2H6_dt = k_silane * p_SiH4 - k_disilane * p_Si2H6;
    dp_H2_dt    = 2 * k_silane * p_SiH4 - k_H2des * p_H2;
    dtheta_H_dt = 0; % update based on adsorption/desorption dynamics
    
    dydt = [dp_SiH4_dt; dp_Si2H6_dt; dp_H2_dt; dtheta_H_dt];
end

