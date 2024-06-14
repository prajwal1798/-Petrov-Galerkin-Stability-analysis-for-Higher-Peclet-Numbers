# Petrov-Galerkin Stability analysis for Higher Peclet Numbers

• Finite Element method for Convection diffusion was formulated using Standard Galerkin, where we observed the generation of weaker diagonal for Advection term. This produced significant instability for Convection dominated problems.

• As a part of treatment for instability, we introduced artificial diffusion by introducing Petrov-Galerkin Method where an additional weight which was indeed a derivative quantity and added to the Advection, which played crucial role in bringing back the stability for advection domination problem.

• Further research was performed on the relaxation function α, where we observed there exists both a critical and an optimal value for the stabilized solution and observed that Optimal value of α performs exactly to Analytical Solution

• Grid Independence was tested for Petrov-Galerkin method and for a velocity of 25m/s, grid independence was attained for nearly 150 Nodes.. 

Swansea University
Zienkiewicz Centre for Computational Engineering
====================================================================
Convection-Diffusion Solver using Finite Element Method (Coursework 2)
=====================================================================

Overview
--------
This MATLAB script provides solutions to the convection-diffusion problem using four different methods:
1. Petrov-Galerkin Method with Optimal Alpha (PG Optimal)
2. Petrov-Galerkin Method with Alpha = 1 (PG Alpha1)
3. Standard Galerkin Method (SG)
4. Comparison of Methods with Analytical Solution (CP)

The script allows the user to select a method, input the number of nodes, and the velocity value. It calculates and displays essential parameters like grid size, wave velocity, Peclet number, number of nodes, and elements. Then, it solves the convection-diffusion equation and plots the solution.

USER GUIDE:

How to Run
--------------------------------------------------------------------------------------------------------------
1. Open MATLAB.
2. Navigate to the directory containing the script.
3. Run the script by typing `convection_diffusion_solver` in the MATLAB command window.
4. Follow the on-screen prompts to select a method and enter the required parameters.

Requirements
--------------------------------------------------------------------------------------------------------------
- MATLAB (The script was designed to be compatible with most recent versions of MATLAB.)

Input Parameters
--------------------------------------------------------------------------------------------------------------
- Number of Nodes: Enter an integer value for the total number of nodes in the domain.
- Velocity Value (m/s): Enter a numerical value for the velocity.

Selection of Method
--------------------------------------------------------------------------------------------------------------
Choose the method to solve the convection-diffusion problem:
1. Enter '1' for the Petrov-Galerkin Method with Optimal Alpha.
2. Enter '2' for the Petrov-Galerkin Method with Alpha = 1.
3. Enter '3' for the Standard Galerkin Method.
4. Enter '4' for a comparison with the Analytical Solution.

The script will then proceed to solve the problem using the selected method and plot the solution.

Note: For method comparison (option 4), the script will plot solutions from both Petrov-Galerkin methods, Standard Galerkin, and the Analytical Solution for direct comparison.

Acknowledgments:
Professor Perumal Nithiarasu for concepts of CFD using FEM, applying stabilization using Petrov-Galerkin methodology & scripting matrix assemblies.

Author:
Prajwal Bharadwaj
2337862
MSc Computational Engineering (2023-25)
Swansea University
Wales, UK


CODE : 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program for the Convection-Diffusion Solver using FEM and applying
% Stabilization techniques 
% 
% Zienkiewicz Centre for Computational Engineering 
% College of Engineering
% Swansea University
%
% MODULE: EGM06- Computational Fluid Dynamics (Assignment - 2)
% Instructor ============= Prof. P Nithiarasu
% Author     ============= Prajwal Bharadwaj (2337862)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function convection_diffusion_solver()
    clc;  
    close all;  

    % UI setup for the users to select:
    disp('Select the method to solve the convection-diffusion problem:');
    disp('1 - Petrov-Galerkin Method with Optimal Alpha (PG Optimal)');
    disp('2 - Petrov-Galerkin Method with Alpha = 1 (PG Alpha1)');
    disp('3 - Standard Galerkin Method (SG)');
    disp('4 - Comparison of methods with Analytical Solution (CP)');
    choice = input('Enter your choice (1, 2, 3, or 4): ');

    % User inputs parameters
    NP = input('Enter the number of nodes: ');
    u = input('Enter the velocity value (m/s): ');
    k = 1; % Diffusion coefficient
    L = 1; % Total domain length
    le = L / (NP - 1); % Length of each element
    pec_nu = (u * le) / (2 * k); % Peclet number

    % Display Parameters
    fprintf('\nGrid Size: %.4fm\n', le);
    fprintf('Wave Velocity: %.1fm/s\n', u);
    fprintf('Peclet Number: %.4f\n', pec_nu);
    fprintf('Number of Nodes: %d\n', NP);
    fprintf('Number of Elements: %d\n\n', NP - 1);
    
    % User choices
    switch choice
        case 1 % Petrov-Galerkin Method with Optimal Alpha
            phi = petrov_galerkin(NP, u, k, L, (u*L/(NP-1))/(2*k)); 
            plot_results(phi, L, 'Petrov-Galerkin Method with Optimal Alpha');

        case 2 % Petrov-Galerkin Method with Alpha = 1
            phi = petrov_galerkin_alpha1(NP, u, k, L); 
            plot_results(phi, L, 'Petrov-Galerkin Method with Alpha = 1');

        case 3 % Standard Galerkin Method
            phi = standard_galerkin(NP, u, k, L);
            plot_results(phi, L, 'Standard Galerkin Method');

        case 4 % Comparison with Analytical Solution
            phi_PG_Optimal = petrov_galerkin(NP, u, k, L, (u*L/(NP-1))/(2*k));
            phi_PG_Alpha1 = petrov_galerkin_alpha1(NP, u, k, L);
            phi_SG = standard_galerkin(NP, u, k, L);
            phi_analytical = analytical_solution(NP, u, k, L);
            compare_methods(phi_PG_Optimal, phi_PG_Alpha1, phi_SG, phi_analytical, L);

        otherwise
            disp('Invalid choice. Exiting.');
            return;
    end
    end

    function plot_results(phi, L, titleText)
    x = linspace(0, L, length(phi));
    plot(x, phi, 'b^-', 'LineWidth', 1, 'MarkerIndices', 1:length(phi));
    grid on;
    xlabel('Domain (x)');
    ylabel('\phi');
    title(titleText);
    end

    function compare_methods(phi_PG_Optimal, phi_PG_Alpha1, phi_SG, phi_analytical, L)
    x = linspace(0, L, length(phi_PG_Optimal));
    figure;   
    plot(x, phi_PG_Optimal, 'r-p', 'LineWidth', 1, 'MarkerSize', 8, 'DisplayName', 'PG Optimal');    
    hold on;    
    plot(x, phi_PG_Alpha1, 'g-^', 'LineWidth', 1, 'MarkerSize', 8, 'DisplayName', 'PG Alpha=1');
    plot(x, phi_SG, 'b-s', 'LineWidth', 1, 'MarkerSize', 8, 'DisplayName', 'Standard Galerkin'); 
    plot(x, phi_analytical, 'k-o', 'LineWidth', 1, 'MarkerSize', 8, 'DisplayName', 'Analytical Solution');    
    hold off;
    grid on;
    legend('Location', 'southwest');
    xlabel('Domain (x)');
    ylabel('\phi (Scalar Quantity)');
    title('Comparison of Methods');
    end

--------------------------------------------------
Petrov_Galerkin:
--------------------------------------------------
    function phi = petrov_galerkin(NP, u, k, L, pec_nu)
    le = L / (NP - 1);                    % Length of each element 
    NE = NP - 1;                          % Number of elements 
    alpha_opt = coth(pec_nu) - 1/pec_nu;  % Optimal Alpha
    
    %% Initialise Global Matrix and Global Forcing function
    gmat = zeros(NP, NP);
    f = zeros(NP, 1);
    
    %% Computing Local and Global Matrix 
    for ie = 1:NE
        emat = zeros(2,2);
        emat(1,1) = -0.5 * u + (k / le) + (alpha_opt * u / 2);
        emat(1,2) = 0.5 * u - (k / le) - (alpha_opt * u / 2);
        emat(2,1) = -0.5 * u - (k / le) - (alpha_opt * u / 2);
        emat(2,2) = 0.5 * u + (k / le) + (alpha_opt * u / 2);

        ip1 = ie;
        ip2 = ie + 1;
        gmat(ip1,ip1) = gmat(ip1,ip1) + emat(1,1);
        gmat(ip1,ip2) = gmat(ip1,ip2) + emat(1,2);
        gmat(ip2,ip1) = gmat(ip2,ip1) + emat(2,1);
        gmat(ip2,ip2) = gmat(ip2,ip2) + emat(2,2);
    end
    
    %% Global Matrix reduction for Applied BC's
    %% phi(1) = 1; phi(0) = 0; 
    % Therefore reduce the first and last row and column repsectively
    gmat_reduced = gmat(2:end-1, 2:end-1);
    f_reduced = f(2:end-1);

    %% Since first column could not be deletd, as the dirichet is not 0, taking it into RHS
    f_reduced(1) = f_reduced(1) - gmat(2, 1) * 1;

    %% Solving reduced solution vector phi
    phi_internal = gmat_reduced \ f_reduced;

    %% Assembling Global solution vector phi
    phi = [1; phi_internal; 0];
    end
-------------------------------------------
Petrov_Galerkin with Optimal Alpha:
-------------------------------------------
    function phi = petrov_galerkin_alpha1(NP, u, k, L)
    le = L / (NP - 1);                    % Length of each element 
    NE = NP - 1;                          % Number of elements 
    alpha_opt = 1;                        % Alpha set to 1
    
    %% Initialise Global Matrix and Global Forcing function
    gmat = zeros(NP, NP);
    f = zeros(NP, 1);
    
    %% Computing Local and Global Matrix 
    for ie = 1:NE
        emat = zeros(2,2);
        emat(1,1) = -0.5 * u + (k / le) + (alpha_opt * u / 2);
        emat(1,2) = 0.5 * u - (k / le) - (alpha_opt * u / 2);
        emat(2,1) = -0.5 * u - (k / le) - (alpha_opt * u / 2);
        emat(2,2) = 0.5 * u + (k / le) + (alpha_opt * u / 2);

        ip1 = ie;
        ip2 = ie + 1;
        gmat(ip1,ip1) = gmat(ip1,ip1) + emat(1,1);
        gmat(ip1,ip2) = gmat(ip1,ip2) + emat(1,2);
        gmat(ip2,ip1) = gmat(ip2,ip1) + emat(2,1);
        gmat(ip2,ip2) = gmat(ip2,ip2) + emat(2,2);
    end
    
    %% Global Matrix reduction for Applied BC's
    gmat_reduced = gmat(2:end-1, 2:end-1);
    f_reduced = f(2:end-1);
    f_reduced(1) = f_reduced(1) - gmat(2, 1) * 1; 

    %% Solving reduced solution vector phi
    phi_internal = gmat_reduced \ f_reduced;

    %% Assembling Global solution vector phi
    phi = [1; phi_internal; 0];
    end

---------------------------------------
Standard_Galerkin:
---------------------------------------

    function phi = standard_galerkin(NP, u, k, L)
    le = L / (NP - 1);            % Length of each element
    NE = NP - 1;                  % Number of elements
    
    %% Initialise Global Matrix and Global Forcing function
    gmat = zeros(NP, NP);
    f = zeros(NP, 1);

    %% Computing Local and Global Matrix 
    for ie = 1:NE
        emat = zeros(2,2);
        emat(1,1) = -0.5 * u + (k / le);
        emat(1,2) = 0.5 * u - (k / le);
        emat(2,1) = -0.5 * u - (k / le);
        emat(2,2) = 0.5 * u + (k / le);

        ip1 = ie;
        ip2 = ie + 1;
        gmat(ip1,ip1) = gmat(ip1,ip1) + emat(1,1);
        gmat(ip1,ip2) = gmat(ip1,ip2) + emat(1,2);
        gmat(ip2,ip1) = gmat(ip2,ip1) + emat(2,1);
        gmat(ip2,ip2) = gmat(ip2,ip2) + emat(2,2);
    end
    

    %% Global Matrix reduction for Applied BC's
    %% phi(1) = 1; phi(0) = 0; 
    % Therefore reduce the first and last row and column repsectively
    gmat_reduced = gmat(2:end-1, 2:end-1);
    f_reduced = f(2:end-1);

    %% Since first column could not be deletd, as the dirichet is not 0, taking it into RHS
    f_reduced(1) = f_reduced(1) - gmat(2, 1) * 1;
    
    %% Solving reduced solution vector phi
    phi_internal = gmat_reduced \ f_reduced;

    %% Assembling Global solution vector phi
    phi = [1; phi_internal; 0];
    end

