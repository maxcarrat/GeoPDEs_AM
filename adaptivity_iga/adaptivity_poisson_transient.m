% ADAPTIVITY_LAPLACE: solve the Poisson transient problem with an adaptive
% isogeometric method based on hierarchical splines and an implicit time
% integration scheme (backward Euler)
%
% [geometry, hmsh, hspace, u, solution_data] = adaptivity_poisson_transient (problem_data, method_data, adaptivity_data, plot_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:            name of the file containing the geometry
%    - nmnn_sides:          sides with Neumann boundary condition (may be empty)
%    - drchlt_sides:        sides with Dirichlet boundary condition
%    - c_diff:              diffusion coefficient (see solve_laplace)
%    - grad_c_diff:         gradient of the diffusion coefficient (if not present, it is taken as zero)
%    - c_cap:               heat capacity
%    - f:                   function handle of the source term
%    - g:                   function for Neumann condition (if nmnn_sides is not empty)
%    - h:                   function for Dirichlet boundary condition (if drchlt_sides is not empty)
%    - time_discretization: time value at each time step
%    - flag_nl:             isNonLinear if true solve using Newton-Raphson
%
%  method_data : a structure with discretization data. It contains the fields:
%    - degree:              degree of the spline functions.
%    - regularity:          continuity of the spline functions.
%    - nsub_coarse:         number of subelements with respect to the geometry mesh (1 leaves the mesh unchanged)
%    - nsub_refine:         number of subelements to be added at each refinement step (2 for dyadic)
%    - nquad:               number of points for Gaussian quadrature rule
%    - space_type:          'simplified' (only children of removed functions) or 'standard' (full hierarchical basis)
%    - truncated:           false (classical basis) or true (truncated basis)
%
%  adaptivity_data: a structure with data for the adaptive method. It contains the fields:
%    - flag:          refinement procedure, based either on 'elements' or on 'functions'
%    - mark_strategy: marking strategy. See 'adaptivity_mark' for details
%    - mark_param:    a parameter to decide how many entities should be marked. See 'adaptivity_mark' for details
%    - max_level:     stopping criterium, maximum number of levels allowed during refinement
%    - max_ndof:      stopping criterium, maximum number of degrees of freedom allowed during refinement
%    - max_nel:       stopping criterium, maximum number of elements allowed during refinement
%    - num_max_iter:  stopping criterium, maximum number of iterations allowed
%    - tol:           stopping criterium, adaptive refinement is stopped when the global error estimator
%                      is lower than tol.
%    - C0_est:        an optional multiplicative constant for scaling the error estimators (default value: 1).
%
%  plot_data: a structure to decide whether to plot things during refinement.
%    - plot_hmesh:        plot the mesh at every iteration
%    - plot_discrete_sol: plot the discrete solution at every iteration
%    - print_info:        display info on the screen on every iteration (number of elements,
%                          number of functions, estimated error, number of marked elements/functions...)
%
% OUTPUT:
%    geometry:      geometry structure (see geo_load)
%    hmsh:          object representing the hierarchical mesh (see hierarchical_mesh)
%    hspace:        object representing the space of hierarchical splines (see hierarchical_space)
%    u:             computed degrees of freedom, at the last iteration.
%    solution_data: a structure with the following fields
%      - iter:       iteration on which the adaptive procedure stopped
%      - ndof:       number of degrees of freedom for each computed iteration
%      - nel:        number of elements for each computed iteration
%      - gest:       global error estimator, for each computed iteration
%      - err_h1s:    error in H1 seminorm for each iteration, if the exact solution is known
%      - err_h1:     error in H1 norm for each iteration, if the exact solution is known
%      - err_l2:     error in L2 norm for each iteration, if the exact solution is known
%      - flag:       a flag with one of the following values:
%          -1: the coefficients for the partition of unity were wrong. This is probably caused by a bug.
%           1: convergence is reached, the global estimator is lower than the given tolerance
%           2: maximum number of iterations reached before convergence.
%           3: maximum number of levels reached before convergence
%           4: maximum number of degrees of freedom reached before convergence
%           5: maximum number of elements reached before convergence
%
%
% For more details about the implementation, see:
%    E. M. Garau, R. Vazquez, Algorithms for the implementation of adaptive
%     isogeometric methods using hierarchical splines, Tech. Report, IMATI-CNR, 2016
%
% For details about the 'simplified' hierarchical space:
%    A. Buffa, E. M. Garau, Refinable spaces and local approximation estimates
%     for hierarchical splines, IMA J. Numer. Anal., (2016)
%
% Copyright (C) 2015, 2016 Eduardo M. Garau, Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.


function  [geometry, hmsh, hspace, u, solution_data] = adaptivity_poisson_transient (problem_data, method_data, adaptivity_data, plot_data)

if (nargin == 3)
    plot_data = struct ('print_info', true, 'plot_hmesh', false, 'plot_discrete_sol', false);
end
if (~isfield (plot_data, 'print_info'))
    plot_data.print_info = true;
end
if (~isfield (plot_data, 'plot_hmesh'))
    plot_data.plot_hmesh = false;
end
if (~isfield (plot_data, 'plot_discrete_sol'))
    plot_data.plot_discrete_sol = false;
end


nel = zeros (1, adaptivity_data.num_max_iter); ndof = nel; gest = nel+NaN;

if (isfield (problem_data, 'graduex'))
    err_h1 = gest;
    err_l2 = gest;
    err_h1s = gest;
end

% Initialization of the hierarchical mesh and space
[hmsh, hspace, geometry] = adaptivity_initialize_laplace (problem_data, method_data);
% Initial solution
number_ts = length(problem_data.time_discretization);
u = ones(hspace.ndof, 1)*problem_data.initial_temperature;
u_last = u;
hspace.dofs = u; 

% ==POST-PROCESSING INITIAL PROBLEM====================================================
if (plot_data.print_info)
    fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Post-Process initial problem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
end
% EXPORT VTK FILE
if (plot_data.print_info); fprintf('\n VTK Post-Process'); end
npts = [plot_data.npoints_x plot_data.npoints_y plot_data.npoints_z];
output_file = sprintf(plot_data.file_name, 0);
sp_to_vtk (hspace.dofs, hspace, geometry, npts(1:hmsh.rdim), output_file, {'solution'}, {'value'})

% Plot in Octave/Matlab
if (plot_data.plot_matlab)
    if (plot_data.print_info); fprintf('\n Octave Post-Process'); end
    if hmsh.ndim == 1
        npts = [plot_data.npoints_x];
        [eu, F] = sp_eval (hspace.dofs, hspace, geometry, npts);
        figure(1000 + 0); plot (squeeze(F(1,:,:)), eu)
    elseif hmsh.ndim == 2
        npts = [plot_data.npoints_x plot_data.npoints_y];
        [eu, F] = sp_eval (hspace.dofs, hspace, geometry, npts);
        figure(1000 + 0); surf (squeeze(F(1,:,:)), squeeze(F(2,:,:)), eu)
    else
        [eu, F] = sp_eval (hspace.dofs, hspace, geometry, npts);
        figure(1000 + 0); surf (squeeze(F(1,:,:)), squeeze(F(2,:,:)), squeeze(F(3,:,:)), eu)
    end
end


%% BACKWARD EULER =========================================================
for itime = 1:number_ts-1
    
    % Initialization of some auxiliary variables
    if ~(isempty(find(plot_data.time_steps_to_post_process==itime, 1)))
        if (plot_data.plot_hmesh)
            fig_mesh = figure(1+itime);
        end
        %     if (plot_data.plot_discrete_sol)
        %         fig_sol = figure(10000+itime);
        %     end
    end
 
    if (plot_data.print_info)
        fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Time step %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n',itime);
    end
    
    
    % ADAPTIVE LOOP
    for iter = 1:adaptivity_data.num_max_iter
        if (plot_data.print_info)
            fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Adaptivity iteration %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n',iter);
        end
        
        if (~hspace_check_partition_of_unity (hspace, hmsh))
            disp('ERROR: The partition-of-the-unity property does not hold.')
            solution_data.flag = -1; break
        end
        
%% SOLVE ==================================================================
        if (plot_data.print_info)
            disp('SOLVE:')
            fprintf('Number of elements: %d. Total DOFs: %d \n', hmsh.nel, hspace.ndof);
        end
        
        if ~(problem_data.flag_nl == 1)
            u = adaptivity_solve_poisson_transient (hmsh, hspace, itime, problem_data);
        else
            [u, problem_data] = adaptivity_solve_poisson_transient_nl (hmsh, hspace, itime, problem_data, plot_data);
        end
        
        nel(iter) = hmsh.nel; ndof(iter) = hspace.ndof;
        
        if (~(isempty(find(plot_data.time_steps_to_post_process==itime, 1))))
            if (plot_data.plot_hmesh)
                fig_mesh = hmsh_plot_cells (hmsh, 20, (fig_mesh) );
            end
        end
        
%% ESTIMATE ===============================================================
        if (plot_data.print_info); fprintf('\n ESTIMATE: \n'); end
        est = adaptivity_estimate_poisson_nl(u, u_last, itime, hmsh, hspace, problem_data, adaptivity_data);
%         est = adaptivity_estimate_gradient(u, itime, hmsh, hspace, problem_data, adaptivity_data);
        gest(iter) = norm (est);
        if (plot_data.print_info); fprintf('Computed error estimator: %f \n', gest(iter)); end
        if (isfield (problem_data, 'graduex'))
            [err_h1(iter), err_l2(iter), err_h1s(iter)] = sp_h1_error (hspace, hmsh, u(:,itime+1), problem_data.uex, problem_data.graduex);
            if (plot_data.print_info); fprintf('Error in H1 seminorm = %g\n', err_h1s(iter)); end
        end
        
% STOPPING CRITERIA -------------------------------------------------------
        if (gest(iter) < adaptivity_data.tol)
            if (plot_data.print_info); disp('Success: The solution converge!!!'); end;
            hspace.dofs = u;
            break;
        elseif (hspace.ndof > adaptivity_data.max_ndof)
            if (plot_data.print_info); disp('Warning: reached the maximum number of DOFs'); end;
            hspace.dofs = u;
            break;
        elseif (hmsh.nel > adaptivity_data.max_nel)
            if (plot_data.print_info); disp('Warning: reached the maximum number of elements'); end;
            hspace.dofs = u; 
            break;
         elseif (iter >= adaptivity_data.num_max_iter)
            if (plot_data.print_info); disp('Warning: reached the maximum number of iterations'); end;
            hspace.dofs = u; 
            break;
        end
%% REFINEMENT =============================================================
% MARK REFINEMENT
        if (plot_data.print_info); disp('MARK REFINEMENT:'); end
        [marked_ref, num_marked_ref] = adaptivity_mark (est, hmsh, hspace, adaptivity_data);
        
        % REFINE
        if (plot_data.print_info)
            fprintf('%d %s marked for refinement \n', num_marked_ref, adaptivity_data.flag);
            disp('REFINE')
        end
        [hmsh, hspace, Cref] = adaptivity_refine (hmsh, hspace, marked_ref, adaptivity_data);
        
        % Project the previous solution mesh onto the next refined mesh
        if (plot_data.print_info); fprintf('\n Project old solution onto refined mesh \n'); end
        
        % project dof onto new mesh
        hspace.dofs = Cref * u;
        u = hspace.dofs;
        % project last time step solution onto new mesh
        u_last = Cref * u_last;
        % project error estimation onto new mesh
        if strcmp(adaptivity_data.flag, 'functions')
            est = Cref * est;
        end;
        
%% COARSENING =============================================================
        % MARK COARSENING
        if (plot_data.print_info); disp('MARK COARSENING:'); end
        [marked_coarse, num_marked_coarse] = marking_for_coarsening (est, hmsh, hspace, adaptivity_data);
        
        if (itime > 1 && ~isempty(marked_coarse))
            % COARSE
            if (plot_data.print_info)
                fprintf('%d %s marked for coarsening \n', num_marked_coarse, adaptivity_data.flag);
                disp('COARSE')
            end
            % Project the previous solution mesh onto the next refined mesh
            if (plot_data.print_info); fprintf('\n Project old solution onto coarsed mesh \n'); end
            % project dofs onto new mesh
            hspace.dofs = u;
            [hmsh_coarse, hspace_coarse, u] = adaptivity_coarsen (hmsh, hspace, marked_coarse, adaptivity_data);
            
            if (plot_data.print_info); fprintf('\n Project last convergent time step \n'); end
            % project last time step solution onto new mesh
            hspace.dofs = u_last;
            [~, ~, u_last] = adaptivity_coarsen (hmsh, hspace, marked_coarse, adaptivity_data);
            
%             if (plot_data.print_info); fprintf('\n Project estimated error \n'); end
%             % project error onto new mesh
%             hspace.dofs = est;
%             [~, ~, est] = adaptivity_coarsen (hmsh, hspace, marked_coarse, adaptivity_data);
            
            hmsh = hmsh_coarse;
            hspace = hspace_coarse;
            hspace.dofs = u;
            
            %             % project error estimation onto new mesh
            %             if strcmp(adaptivity_data.flag, 'functions')
            %                 est = Ccoarse * est;
            %             end;
        end
        
%% ========================================================================
        
        if (plot_data.adaptivity)
            % ==POST-PROCESSING====================================================
            if (plot_data.print_info)
                fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Post-Process adaptivity = %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n', iter);
            end
            % EXPORT VTK FILE
            if (plot_data.print_info); fprintf('\n VTK Post-Process'); end
            npts = [plot_data.npoints_x plot_data.npoints_y plot_data.npoints_z];
            output_file = sprintf(plot_data.file_name, iter);
            sp_to_vtk (hspace.dofs, hspace, geometry, npts(1:hmsh.rdim), output_file, {'solution', 'gradient'}, {'value', 'gradient'})
            if strcmp(adaptivity_data.flag, 'functions')
                output_file_est = sprintf(plot_data.file_name_err, iter);
                sp_to_vtk (est, hspace, geometry, npts(1:hmsh.rdim), output_file_est, {'error'})
            end
            
            % Plot in Octave/Matlab
            if (plot_data.plot_matlab)
                if (plot_data.print_info); fprintf('\n Octave Post-Process'); end
                [eu, F] = sp_eval (hspace.dofs, hspace, geometry, npts);
                if numel(hmsh.rdim) == 1
                    figure(1000 + 0); plot (squeeze(F(1,:,:)), eu)
                elseif numel(hmsh.rdim) == 2
                    figure(1000 + 0); surf (squeeze(F(1,:,:)), squeeze(F(2,:,:)), eu)
                else
                    figure(1000 + 0); surf (squeeze(F(1,:,:)), squeeze(F(2,:,:)), squeeze(F(3,:,:)), eu)
                end
            end
        end
        
    end % END ADAPTIVITY LOOP
    
    %update last time step solution
    u_last = u;
    
    solution_data.iter = iter;
    solution_data.gest = gest(1:iter);
    solution_data.ndof = ndof(1:iter);
    solution_data.nel  = nel(1:iter);
    if (exist ('err_h1s', 'var'))
        solution_data.err_h1s = err_h1s(1:iter);
        solution_data.err_h1 = err_h1(1:iter);
        solution_data.err_l2 = err_l2(1:iter);
    end
    
    if ~(isempty(find(plot_data.time_steps_to_post_process==itime, 1)))
        % ==POST-PROCESSING====================================================
        if (plot_data.print_info)
            fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Post-Process time step = %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n', itime);
        end
        % EXPORT VTK FILE
        if (plot_data.print_info); fprintf('\n VTK Post-Process'); end
        npts = [plot_data.npoints_x plot_data.npoints_y plot_data.npoints_z];
        output_file = sprintf(plot_data.file_name, itime);
        sp_to_vtk (hspace.dofs, hspace, geometry, npts(1:hmsh.rdim), output_file, {'solution', 'gradient'}, {'value', 'gradient'})
        if strcmp(adaptivity_data.flag, 'functions')
            output_file_est = sprintf(plot_data.file_name_err, itime);
            sp_to_vtk (est, hspace, geometry, npts(1:hmsh.rdim), output_file_est, {'error'})
        end
        
        % Plot in Octave/Matlab
        if (plot_data.plot_matlab)
            if (plot_data.print_info); fprintf('\n Octave Post-Process'); end
            if hmsh.ndim == 1
                npts = [plot_data.npoints_x];
                [eu, F] = sp_eval (u, hspace, geometry, npts);
                figure(1000 + itime); plot (squeeze(F(1,:,:)), eu)
            elseif hmsh.ndim == 2
                npts = [plot_data.npoints_x plot_data.npoints_y];
                [eu, F] = sp_eval (u, hspace, geometry, npts);
                figure(1000 + itime); surf (squeeze(F(1,:,:)), squeeze(F(2,:,:)), eu)
            else
                [eu, F] = sp_eval (u, hspace, geometry, npts);
                figure(1000 + itime); surf (squeeze(F(1,:,:)), squeeze(F(2,:,:)), squeeze(F(3,:,:)), eu)
            end
        end
    end
end % END BACKWARD EULER LOOP

end