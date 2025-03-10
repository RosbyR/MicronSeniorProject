function polysilicon_growth()
    % Retrieve constants and boundary layer conditions
    [R, T, E_decom, A_decom, SiH4_conc, tspan, P, delta, D_SiH4, D_SiH2, D_H2, ...
        A_SiH4, E_SiH4, A_SiH2, E_SiH2, A_H2des, E_H2des, A_SiH2_growth, E_SiH2_growth, Rad, A, Thick, DS, MMs] = Micron_constants1();
    
    [SiH4_b, SiH2_b, H2_b] = Micron_gasphase();
    
    % Ensure values are scalars
    SiH4_b = SiH4_b(1);
    SiH2_b = SiH2_b(1);
    H2_b   = H2_b(1);
    
    % Surface reaction rate constants
    k_SiH2_s = A_SiH2_growth  * exp(E_SiH2_growth  / (R*T));
    
    % Initial surface silicon concentration
    Si_surface = 0; % Initial condition (arbitrary units)

    % Solve ODE for silicon growth
    options = odeset('RelTol',1e-6,'AbsTol',1e-9);
    [t, Si_growth] = ode15s(@(t, Si_surface) growth_ODE(t, Si_surface, k_SiH2_s, SiH4_b, SiH2_b), tspan, Si_surface, options);
    
    Film_thick = (Si_growth/MMs) * (DS) * (A) * (1/1e7) * 2 * 100
    
    % Display results
   
    
    % Plot growth over time
    figure;
  
    plot(t, Si_growth, 'b-', 'LineWidth', 2);
  
    xlabel('Time (s)');
    ylabel('Polysilicon Growth mol');
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
    % Deposition contribution from silane and disilane
    Si_deposit = k_SiH2_s * SiH2_b;
    
    % Total silicon growth rate
    dSi_dt = Si_deposit;
end
