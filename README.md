# Micron Senior Design Project: Modeling and Optimization of a CVD Chamber for Semiconductor Manufacturing
This project was conducted in collaboration with Micron Technologies to design and optimize a Chemical Vapor Deposition (CVD) chamber for semiconductor manufacturing. Using MATLAB-based modeling, we extracted and optimized key design variables such as pressure, temperature, precursor composition, and deposition time. The model was validated against experimental data provided by Micron, and then applied to guide design decisions for improved chamber performance.
---

## Repository Structure  
- **Micron Project.pdf** - Full copy of the teams design report submitted to Micron
- **Micron_constants1.m** – Defines constants and physical parameters used across the model  
- **Micron_deposition.m** – Core deposition kinetics and simulation  
- **Micron_exhaust.m** – Exhaust modeling and mass balance considerations  
- **Micron_minimize.m** – Optimization of key process variables  
- **plot_minimize.m** – Visualization of optimization results  
- **Micron_tornado.m** – Tornado plots for sensitivity analysis  
- **micron_bigplot.m** – Combined visualization of parameter sweeps  
- **README.md** – This document  

---

## Methodology  
1. Developed a MATLAB-based kinetic model of deposition behavior.  
2. Validated simulation outputs against Micron’s provided experimental data.  
3. Performed parameter sweeps across temperature, pressure, and gas composition.  
4. Applied optimization routines to identify conditions that maximized uniformity and minimized defect risk.  
5. Conducted sensitivity analysis to quantify the effect of each parameter.  

---

## Key Results  
- Verified deposition rates matched experimental data within 2%.  
- Identified optimal chamber conditions with a wafer uniformity of 3% range/mean.  
- Delivered a full design study and final recommendations to Micron engineers.  

---

## How to Run  
From MATLAB:  

```matlab
% Run the core deposition model
Micron_deposition  

% Perform optimization
Micron_minimize  

% Generate tornado sensitivity plots
Micron_tornado  
