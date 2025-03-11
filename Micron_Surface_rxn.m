function [SiH4_s, SiH2_s, H2_s] = Micron_Surface_rxn()
    % Retrieve constants in correct order:
    [R, T, E_decom, A_decom, SiH4_conc, H2_conc, tspan, P, delta, D_SiH4, D_SiH2, D_H2, ...
     A_SiH4, E_SiH4, A_SiH2, E_SiH2, A_H2des, E_H2des] = Micron_constants1();
     
    % Get boundary layer concentrations from Micron_gasphase
    [SiH4_b, SiH2_b, H2_b] = Micron_gasphase();
    % Ensure scalar values:
    SiH4_b = SiH4_b(1);
    SiH2_b = SiH2_b(1);
    H2_b   = H2_b(1);
    
    % Compute rate constants based on first-order kinetics
    k_silane   = A_SiH4 * exp(-E_SiH4 / (R*T));
    k_disilane = A_SiH2 * exp(-E_SiH2 / (R*T));
    k_H2des    = A_H2des * exp(-E_H2des / (R*T));
    
    % Compute effective diffusivities for flux
    k_diff_SiH4 = D_SiH4 / delta;
    k_diff_SiH2 = D_SiH2 / delta;
    k_diff_H2   = D_H2   / delta;
    
    % Set initial conditions for the ODE
    % State vector: [SiH4_b; SiH4_s; SiH2_b; SiH2_s; H2_b; H2_s]
    y0 = [SiH4_b; 0; SiH2_b; 0; H2_b; 0];
    
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
    [t, y] = ode15s(@(t, y) SurfaceRxnODE(t, y, k_diff_SiH4, k_silane, ...
                          k_diff_SiH2, k_disilane, k_diff_H2, k_H2des), tspan, y0, options);
    
    % Return the final surface concentrations:
    SiH4_s = y(end, 2);
    SiH2_s = y(end, 4);
    H2_s   = y(end, 6);
    
    % Optional: display final surface concentrations for debugging
    fprintf('Final surface concentrations:\n');
    fprintf('  SiH4: %g\n', SiH4_s);
    fprintf('  SiH2: %g\n', SiH2_s);
    fprintf('  H2:   %g\n', H2_s);
end

function dydt = SurfaceRxnODE(~, y, k_diff_SiH4, k_silane, k_diff_SiH2, k_disilane, k_diff_H2, k_H2des)
    % Unpack state variables:
    SiH4_b = y(1);
    SiH4_s = y(2);
    SiH2_b = y(3);
    SiH2_s = y(4);
    H2_b   = y(5);
    H2_s   = y(6);
    
    % Compute diffusive fluxes for each species:
    flux_SiH4 = k_diff_SiH4 * (SiH4_b - SiH4_s);
    flux_SiH2 = k_diff_SiH2 * (SiH2_b - SiH2_s);
    flux_H2   = k_diff_H2   * (H2_b   - H2_s);
    
    % Total flux (used for normalization):
    flux_tot = flux_SiH4 + flux_SiH2 + flux_H2;
    %Inhibition factor
    beta = 3;  % Adjust as needed
    % Boundary layer equations (gas-phase changes):
    dSiH4_b_dt = -flux_SiH4;
    dSiH2_b_dt = -flux_SiH2;
    dH2_b_dt   = -flux_H2;
    
    % Surface equations (absorbed species):
    dSiH4_s_dt = (flux_SiH4 - k_silane * SiH4_s) / flux_tot;
    dSiH2_s_dt = (flux_SiH2 - k_disilane * SiH2_s) / (flux_tot + beta * H2_s);
    dH2_s_dt   = (flux_H2   - k_H2des   * H2_s) / flux_tot;
    
    % Pack derivatives into a column vector:
    dydt = [dSiH4_b_dt; dSiH4_s_dt; dSiH2_b_dt; dSiH2_s_dt; dH2_b_dt; dH2_s_dt];
end
