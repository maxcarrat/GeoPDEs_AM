% PHYSICAL DATA OF THE PROBLEM
clear problem_data
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_square.txt';
% Generate Output folder
problem_output.folder = 'output';
mkdir(problem_output.folder);

% Set non-linear flag
% if true the non-linear solver is used
problem_data.flag_nl = true;

% Non-linear analysis
problem_data.Newton_tol = 1.0e-06;
problem_data.num_Newton_iter = 5;
problem_data.non_linear_convergence_flag = 0;

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [];

% Physical parameters
problem_data.c_diff  = @conductivity;
problem_data.grad_c_diff = @conductivity_der;
problem_data.c_cap = @capacity;

% Time discretization
n_time_steps = 10;
problem_data.time_discretization = linspace(0.0, 100.0, n_time_steps + 1);

% Heat Source path
x_begin = 0.25;
x_end = 0.75;
x_path = linspace(x_begin, x_end, n_time_steps+1);
y_begin = 0.25;
y_end = 0.75;
y_path = linspace(y_begin, y_end, n_time_steps+1);
problem_data.path = [x_path', y_path'];

% Source and boundary terms
problem_data.f = moving_heat_source;        % Body Load
problem_data.h = [];                        % Dirichlet Boundaries

% CHOICE OF THE DISCRETIZATION PARAMETERS (Coarse mesh)
clear method_data
method_data.degree      = [2 2];        % Degree of the splines
method_data.regularity  = [1 1];        % Regularity of the splines
method_data.nsub_coarse = [5 5];        % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nsub_refine = [2 2];        % Number of subdivisions for each refinement
method_data.nquad       = [3 3];        % Points for the Gaussian quadrature rule
method_data.space_type  = 'standard';   % 'simplified' (only children functions) or 'standard' (full basis)
method_data.truncated   = 0;            % 0: False, 1: True

% ADAPTIVITY PARAMETERS
clear adaptivity_data
adaptivity_data.flag = 'elements';
% adaptivity_data.flag = 'functions';
adaptivity_data.C0_est = 1.0;
adaptivity_data.mark_param = 0.75;
adaptivity_data.mark_strategy = 'GRAD';
adaptivity_data.max_level = 2;
adaptivity_data.max_ndof = 15000;
adaptivity_data.num_max_iter = 5;
adaptivity_data.max_nel = 15000;
adaptivity_data.tol = 1e-03;

% GRAPHICS
plot_data.plot_hmesh = true;
plot_data.plot_discrete_sol = false;
plot_data.print_info = true;
plot_data.time_steps_to_post_process = [1,5,10];%linspace(1,10);%
plot_data.file_name = strcat(problem_output.folder, '/poisson_adaptivity_travelling_heat_source_error2_%d.vts');

[geometry, hmsh, hspace, u, solution_data] = adaptivity_poisson_transient(problem_data, method_data, adaptivity_data, plot_data);