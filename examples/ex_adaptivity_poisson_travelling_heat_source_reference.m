% PHYSICAL DATA OF THE PROBLEM
close all
clear problem_data
clc
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_2DPrintingLayer.txt';
% Generate Output folder
problem_output.folder = '~/workspace/GeoPDEs_AM.git/trunk/output';
mkdir(problem_output.folder);

% Set non-linear flag
% if true the non-linear solver is used
problem_data.flag_nl = false;
problem_data.lumped = true;
problem_data.non_linear_convergence_flag= true;

% Non-linear analysis
problem_data.Newton_tol = 1.0e-05;
problem_data.num_Newton_iter = 20;

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [];

% Physical parameters
problem_data.c_diff  = @(x,y) ones(size(x)) * 29.0; % 
problem_data.grad_c_diff = @(x,y) cat (1,reshape (zeros(size(x)), [1, size(x)]), reshape (zeros(size(x)), [1, size(x)])); % 
problem_data.c_cap =  @(x,y) ones(size(x)) * 7820 * 600; %  
problem_data.initial_temperature = 20;                  %[°C]


% Time discretization
n_time_steps = 20;
time_end = 20;
problem_data.time_discretization = linspace(0.0, time_end, n_time_steps + 1);

% Heat Source path
x_begin = 0.05;
x_end = 0.15;
x_path = linspace(x_begin, x_end, n_time_steps+1);
y_begin = 0.05;
y_end = 0.05;
y_path = linspace(y_begin, y_end, n_time_steps+1);
problem_data.path = [x_path', y_path'];

% Source and boundary terms
problem_data.f = @(x,y,path_x,path_y)  3.0 * 3.0 * 583000.0/(pi * 0.015 * 0.010 ) * ...
    exp(- 3*(x-path_x).^2/0.015.^2 - 3*(y-path_y).^2/0.010.^2 );       % Body Load
problem_data.h = [];                        % Dirichlet Boundaries
problem_data.g = [];                        % Neumann Boundaries

% CHOICE OF THE DISCRETIZATION PARAMETERS (Coarse mesh)
clear method_data
method_data.degree      = [2 2];        % Degree of the splines
method_data.regularity  = [1 1];        % Regularity of the splines
method_data.nsub_coarse = [2^5 2^5];        % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nsub_refine = [1 1];        % Number of subdivisions for each refinement
method_data.nquad       = [3 3];        % Points for the Gaussian quadrature rule
method_data.space_type  = 'standard';   % 'simplified' (only children functions) or 'standard' (full basis)
method_data.truncated   = 1;            % 0: False, 1: True

% ADAPTIVITY PARAMETERS
clear adaptivity_data
adaptivity_data.flag = 'elements';
adaptivity_data.C0_est = 1.0;
adaptivity_data.mark_param = 0.75;
adaptivity_data.mark_param_coarsening = 0.1;
adaptivity_data.mark_neighbours = false;
adaptivity_data.crp = 1.5;                     %coarsening relaxation parameter
adaptivity_data.mark_strategy = 'MS';
adaptivity_data.radius = [0.015, 0.01];
adaptivity_data.max_level = 5;
adaptivity_data.max_ndof = 15000;
adaptivity_data.num_max_iter = 1;
adaptivity_data.max_nel = 15000;
adaptivity_data.tol = 2.0e-01;

% GRAPHICS
plot_data.plot_hmesh = false;
plot_data.adaptivity = false;
plot_data.plot_discrete_sol = false;
plot_data.print_info = true;
plot_data.plot_matlab = false;
plot_data.time_steps_to_post_process = [1, 5, 10, 15, 20]; %linspace(1,n_time_steps,n_time_steps); % [1, 10, 20, 30, 40, 50, 100, 150, 200, 250, 300, 350, 400]; %      
plot_data.file_name = strcat(problem_output.folder, '/REFERENCE_poisson_adaptivity_Fachinotti20_travelling_heat_source_lumped_2D_%d.vts');
plot_data.file_name_mesh = strcat(problem_output.folder, '/REFERENCE_poisson_adaptivity_Fachinotti20_travelling_heat_source_lumped_2D_mesh_%d');
plot_data.npoints_x = 201;        %number of points x-direction in post-processing
plot_data.npoints_y = 101;        %number of points y-direction in post-processing
plot_data.npoints_z = 1;          %number of points z-direction in post-processing

[geometry, hmsh, hspace, u, solution_data] = adaptivity_poisson_transient(problem_data, method_data, adaptivity_data, plot_data);