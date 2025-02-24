[R, T, k1, SiH4_conc, tspan, P, delta, D_SiH4, D_SiH2, D_H2, ...
 A_SiH4, E_SiH4, A_SiH2, E_SiH2, A_H2des, E_H2des] = Micron_constants1();


[SiH4_b, SiH2_b, H2_b] = Micron_gasphase();
% Ensure that the outputs are scalars:
SiH4_b = SiH4_b(1);
SiH2_b = SiH2_b(1);
H2_b   = H2_b(1);



% Compute rate constants based on first order kinetics
k_silane   = A_SiH4  * exp(-E_SiH4  / (R*T));
k_disilane = A_SiH2  * exp(-E_SiH2  / (R*T));
k_H2des    = A_H2des * exp(-E_H2des / (R*T));

% Initial conditions (column vector)
y0 = [SiH4_b; SiH2_b; H2_b];

% Use the corrected parameter order
[t, y] = ode15s(@(t, y) SurfaceRxnODE(t, y, R, k_silane, k_disilane, k_H2des), tspan, y0);

% Display final concentrations
disp('Final surface concentrations:');
disp(['SiH4_s: ', num2str(y(end, 1))]);
disp(['SiH2_s: ', num2str(y(end, 2))]);
disp(['H2_s:   ', num2str(y(end, 3))]);

% ODE function definition with debugging
function dydt = SurfaceRxnODE(~, y, R, k_silane, k_disilane, k_H2des)
    % Unpack current concentrations
    p_SiH4 = y(1);
    p_SiH2 = y(2);
    p_H2   = y(3);
    
    % Compute derivatives (using scalar multiplication)
    dp_SiH4_dt = -k_silane * p_SiH4;  
    dp_SiH2_dt = -k_disilane * p_SiH2;
    dp_H2_dt   = -k_H2des  * p_H2;
    
    % Concatenate into a column vector
    dydt = [dp_SiH4_dt; dp_SiH2_dt; dp_H2_dt];
end
