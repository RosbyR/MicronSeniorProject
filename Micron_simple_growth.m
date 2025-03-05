function polysilicon_growth()
    % Retrieve constants and boundary layer conditions
    [R, T, k1, SiH4_conc, tspan, P, delta, D_SiH4, D_SiH2, D_H2, ...
        A_SiH4, E_SiH4, A_SiH2, E_SiH2, A_H2des, E_H2des] = Micron_constants1();
    
    [SiH4_b, SiH2_b, H2_b] = Micron_gasphase();
    
    % Ensure values are scalars
    SiH4_b = SiH4_b(1);
    SiH2_b = SiH2_b(1);
    H2_b   = H2_b(1);
    
    % Surface reaction rate constants
    k_silane   = A_SiH4  * exp(-E_SiH4  / (R*T));
    k_disilane = A_SiH2  * exp(-E_SiH2  / (R*T));
    
    % Initial surface silicon concentration
    Si_surface = 0; % Initial condition (arbitrary units)
    
    % Solve ODE for silicon growth
    options = odeset('RelTol',1e-6,'AbsTol',1e-9);
    [t, Si_growth] = ode15s(@(t, Si_surface) growth_ODE(t, Si_surface, k_silane, k_disilane, SiH4_b, SiH2_b), tspan, Si_surface, options);
    
    % Display results
   
    
    % Plot growth over time
    figure;
    plot(t, Si_growth, 'b-', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Polysilicon Growth mol/m^3');
    title('Polysilicon Growth over Time');
    grid on;
end

function dSi_dt = growth_ODE(~, Si_surface, k_silane, k_disilane, SiH4_b, SiH2_b)
    % Deposition contribution from silane and disilane
    SiH4_deposit = k_silane * SiH4_b;
    SiH2_deposit = k_disilane * SiH2_b;
    
    % Total silicon growth rate
    dSi_dt = SiH4_deposit + SiH2_deposit;
end
