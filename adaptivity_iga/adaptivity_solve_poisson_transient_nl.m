% ADAPTIVITY_SOLVE_LAPLACE: assemble and solve the linear system for
% non-linear Poisson problem, using hierarchical spaces and implicit
% time integartion scheme (backward Euler)
%
% The function solves the non-linear problem
%
%  c(x,u) d(u)/dt  - div ( epsilon(x,u) grad (u)) = f(t)    in Omega = F((0,1)^n)
%                              epsilon(x,u) du/dn = g       on Gamma_N
%                                               u = h       on Gamma_D
%
% USAGE:
%
% u = adaptivity_solve_poisson_transient (hmsh, hspace, method_data)
%
% INPUT:
%
%   hmsh:   object representing the hierarchical mesh (see hierarchical_mesh)
%   hspace: object representing the space of hierarchical splines (see hierarchical_space)
%   u_prev: solution vector of the last iteration
%   problem_data: a structure with data of the problem. For this function, it must contain the fields:
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - c_diff:       temperature dependent diffusion coefficient (epsilon in the equation)
%    - c_cap:        temperature dependent heat capacity (c in the equation)
%    - f:            function handle of the source term
%    - g:            function for Neumann condition (if nmnn_sides is not empty)
%    - h:            function for Dirichlet boundary condition
%
% OUTPUT:
%
%   u: computed degrees of freedom
%
% Copyright (C) 2015 Eduardo M. Garau, Rafael Vazquez
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

function [u, problem_data] = adaptivity_solve_poisson_transient_nl (hmsh, hspace, time_step, problem_data, plot_data)

iter = 0;

% Set initial solution as the last convergent solution
if ~isempty(hspace.dofs)
    u_0 = hspace.dofs;
else
    u_0 = zeros(hspace.ndof, 1);
end

u = u_0;
res = 1.0;

while (norm(res) > problem_data.Newton_tol)
    iter = iter + 1;
    
    if (plot_data.print_info)
        fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Newton-Raphson iteration %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n',iter);
    end
    
    stiff_mat = op_gradu_gradv_nl_hier (hspace, hspace, hmsh, problem_data.c_diff, u);
    mass_mat = op_u_v_nl_hier(hspace, hspace, hmsh, problem_data.c_cap, u, u_0);
    %lumped mass matrix
    if problem_data.lumped
        mass_mat = diag(sum(mass_mat));
    end
    rhs = op_f_v_time_hier (hspace, hmsh, problem_data, time_step);
    
    % Apply Neumann boundary conditions
    if (~isfield (struct (hmsh), 'npatch')) % Single patch case
        for iside = problem_data.nmnn_sides
            if (hmsh.ndim > 1)
                % Restrict the function handle to the specified side, in any dimension, gside = @(x,y) g(x,y,iside)
                gside = @(varargin) problem_data.g(varargin{:},iside);
                dofs = hspace.boundary(iside).dofs;
                rhs(dofs) = rhs(dofs) + op_f_v_hier (hspace.boundary(iside), hmsh.boundary(iside), gside);
            else
                if (iside == 1)
                    x = hmsh.mesh_of_level(1).breaks{1}(1);
                else
                    x = hmsh.mesh_of_level(1).breaks{1}(end);
                end
                sp_side = hspace.boundary(iside);
                rhs(sp_side.dofs) = rhs(sp_side.dofs) + problem_data.g(x,iside);
            end
        end
    else % Multipatch case
        boundaries = hmsh.mesh_of_level(1).boundaries;
        Nbnd = cumsum ([0, boundaries.nsides]);
        for iref = problem_data.nmnn_sides
            iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
            gref = @(varargin) problem_data.g(varargin{:},iref);
            rhs_nmnn = op_f_v_hier (hspace.boundary, hmsh.boundary, gref, iref_patch_list);
            rhs(hspace.boundary.dofs) = rhs(hspace.boundary.dofs) + rhs_nmnn;
        end
    end
    
    % Apply Dirichlet boundary conditions
    if ~isempty(problem_data.h)
        u = zeros (hspace.ndof, 1);
        [u_dirichlet, dirichlet_dofs] = sp_drchlt_l2_proj (hspace, hmsh, problem_data.h, problem_data.drchlt_sides);
        u(dirichlet_dofs) = u_dirichlet;
        
        int_dofs = setdiff (1:hspace.ndof, dirichlet_dofs);
        rhs(int_dofs) = rhs(int_dofs) - stiff_mat(int_dofs, dirichlet_dofs)*u(dirichlet_dofs) + mass_mat(int_dofs, int_dofs)*u_prev(int_dofs);
        
        % Solve the linear system
        delta_t = problem_data.time_discretization(time_step)-problem_data.time_discretization(time_step+1);
        lhs = mass_mat(int_dofs, int_dofs) - stiff_mat(int_dofs, int_dofs) * delta_t;
        u(int_dofs) =  lhs\ rhs(int_dofs);
    else
        
        delta_t = problem_data.time_discretization(time_step+1) - problem_data.time_discretization(time_step);
        
        %residuum
        if (plot_data.print_info)
            fprintf('\n \t Assembly residuum vector');
        end
        res = rhs * delta_t - stiff_mat * u * delta_t +...
            mass_mat * u_0 - mass_mat * u;
        
        %jacobian of the residuum
        if (plot_data.print_info)
            fprintf('\n \t Assembly tangent matrix');
        end
        J = stiff_mat * delta_t + mass_mat;
        
        % compute temperature update
        if (plot_data.print_info)
            fprintf('\n \t Evaluate solution increment');
        end
        
        delta_u = J \ res;
        u = u + delta_u;
        
        if (plot_data.print_info)
            fprintf('\n res norm = %d \n',norm(res));
        end
        
    end
    
    
    % STOPPING CRITERIA
    if (iter == problem_data.num_Newton_iter)
        disp('Warning: Reached max number of iterations')
        break;
    end
    
end % END NEWTON-RAPHSON ITERATION LOOP

% UPDATE ADAPTIVITY NON-LINEAR FLAG
if (norm(res) <= problem_data.Newton_tol)
    problem_data.non_linear_convergence_flag = 1;
end

end