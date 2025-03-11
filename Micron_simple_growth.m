function polysilicon_growth()
    % Retrieve constants and boundary layer conditions
    [R, T, E_decom, A_decom, SiH4_conc, H2_conc, tspan, P, delta, D_SiH4, D_SiH2, D_H2, ...
        A_SiH4, E_SiH4, A_SiH2, E_SiH2, A_H2des, E_H2des, A_SiH2_growth, E_SiH2_growth, ...
        Rad, A, Thick, DS, MMs, Ratio] = Micron_constants1();
    
    % Get surface concentrations from Micron_Surface_rxn
    [SiH4_s, SiH2_s, H2_s] = Micron_Surface_rxn();
    
    % Ensure values are scalars
    SiH4_s = SiH4_s(1);
    SiH2_s = SiH2_s(1);
    H2_s   = H2_s(1);
    
    % Calculate the surface reaction rate constant
    k_SiH2_s = A_SiH2_growth * exp(E_SiH2_growth / (R*T));
    
    % Display diagnostic values to ensure they are as expected
    fprintf('Diagnostic values at t = 0:\n');
    fprintf('SiH4_s: %g\n', SiH4_s);
    fprintf('SiH2_s: %g\n', SiH2_s);
    fprintf('H2_s:   %g\n', H2_s);
    fprintf('k_SiH2_s: %g\n', k_SiH2_s);
    
    % Initial surface silicon concentration (state variable)
    Si_surface = .0001; % (arbitrary units)
    y0 = Si_surface;
    
    % Evaluate the derivative at t = 0 to check for any anomalies
    dSi_dt0 = growth_ODE(0, Si_surface, k_SiH2_s, SiH4_s, SiH2_s);
    fprintf('Initial derivative dSi_dt: %g\n', dSi_dt0);
    
    % Solve ODE for silicon growth
    options = odeset('RelTol',1e-6,'AbsTol',1e-9);
    [t, Si_growth] = ode15s(@(t, Si_surface) growth_ODE(t, Si_surface, k_SiH2_s, SiH4_s, SiH2_s), tspan, y0, options);
    
    Film_thick = ((Si_growth/MMs) * (DS) * (A) * (1/1e7) * 2 * 100)/2;
    fprintf('Ratio of H2 to Silane in Precursor Gas %g\n', Ratio);
    % Plot growth over time
    figure;
    plot(t, Si_growth, 'b-', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Polysilicon Growth (mol)');
    title('Polysilicon Growth over Time');
    grid on;
    
    figure;
    plot(t, Film_thick, 'b-', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Polysilicon Film Thickness (Angstrom)');
    title('Polysilicon Film Growth over Time');
    grid on;
end

function dSi_dt = growth_ODE(~, Si_surface, k_SiH2_s, SiH4_b, SiH2_b)
    % Compute the deposition contribution from disilane
    Si_deposit = k_SiH2_s * SiH2_b;
    
    % Total silicon growth rate (assumed to be solely based on Si_deposit)
    dSi_dt = Si_deposit;
end
