% PHYSICAL DATA OF THE PROBLEM
clear problem_data
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_2DPrintingLayer.txt';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [];

% Physical parameters
problem_data.c_diff  = @(x, y) ones(size(x)) * 29.0;
problem_data.grad_c_diff = @(x, y) cat (1, ...
    reshape (zeros(size(x)), [1, size(x)]), ...
    reshape (zeros(size(x)), [1, size(x)]));
problem_data.c_cap = @(x, y) ones(size(x)) * 4.629e+06;

% Time discretization
n_time_steps = 20;
problem_data.time_discretization = linspace(0.0, 20.0, n_time_steps + 1);

% Heat Source path
x_begin = 0.05;
x_end = 0.15;
x_path = linspace(x_begin, x_end, n_time_steps+1);
y_begin = 0.05;
y_end = 0.05;
y_path = linspace(y_begin, y_end, n_time_steps+1);
problem_data.path = [x_path', y_path'];

% Source and boundary terms
problem_data.moving_heat_source = moving_thermal_heat_source;
problem_data.f = problem_data.moving_heat_source;
problem_data.h = [];
problem_data.initial_temperature = 20;                  %[°C]

% CHOICE OF THE DISCRETIZATION PARAMETERS (Coarse mesh)
clear method_data
method_data.degree     = [2 2];       % Degree of the splines
method_data.regularity = [1 1];       % Regularity of the splines
method_data.nsub       = [20 20];     % Number of subdivisions
method_data.nquad      = [3 3];       % Points for the Gaussian quadrature rule

% GRAPHICS
plot_data.plot_hmesh = true;
plot_data.plot_discrete_sol = false;
plot_data.print_info = true;

[geometry, msh, space, u] = solve_poisson_transient (problem_data, method_data);

% 4) POST-PROCESSING
vtk_pts = {linspace(0, 1,100), linspace(0, 1, 100)};

% 4.1) PLOT IN MATLAB. COMPARISON WITH THE EXACT SOLUTION

[eu, F] = sp_eval (u, space, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
surf (X, Y, eu)
title ('Numerical solution'), axis tight

% 4.2) EXPORT TO PARAVIEW

output_file = 'output/Transient_Poisson_moving_source';

fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

