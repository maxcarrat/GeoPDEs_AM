% SOLVE_POISSON_TRANSIENT: assemble and solve the linear system for
% Poisson problem, using tensor product spaces and implicit time integration
% scheme (backward Euler)
%
% The function solves the diffusion problem
%
%  c d(u)/dt  - div ( epsilon(x) grad (u)) = f(t)    in Omega = F((0,1)^n)
%                             epsilon(x) du/dn = g       on Gamma_N
%                                            u = h       on Gamma_D
%
% USAGE:
%
% u = solve_poisson_transient (hmsh, hspace, method_data)
%
% INPUT:
%
%   problem_data: a structure with data of the problem. For this function, it must contain the fields:
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - c_diff:       diffusion coefficient (epsilon in the equation)
%    - c_cap:        heat capacity (c in the equation)
%    - f:            function handle of the source term
%    - g:            function for Neumann condition (if nmnn_sides is not empty)
%    - h:            function for Dirichlet boundary condition
%
%  method_data : a structure with discretization data. Its fields are:
%    - degree:     degree of the spline functions.
%    - regularity: continuity of the spline functions.
%    - nsub:       number of subelements with respect to the geometry mesh
%                   (nsub=1 leaves the mesh unchanged)
%    - nquad:      number of points for Gaussian quadrature rule
%
% OUTPUT:
%
%  geometry: geometry structure (see geo_load)
%  msh:      mesh object that defines the quadrature rule (see msh_cartesian)
%  space:    space object that defines the discrete space (see sp_scalar)
%  u:        the computed degrees of freedom
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

function [geometry, msh, space, u] = solve_poisson_transient (problem_data, method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
    eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
    eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

% Construct geometry structure
geometry  = geo_load (geo_name);

[knots, zeta] = kntrefine (geometry.nurbs.knots, nsub-1, degree, regularity);

% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_cartesian (zeta, qn, qw, geometry);

% Construct space structure
space    = sp_bspline (knots, degree, msh);

% Initial solution
u_prev = zeros(space.ndof, 1);

% BACKWARD EULER LOOP
for itime = 1:length(problem_data.time_discretization)-1
    fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Time step %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n',itime);
   
    stiff_mat = op_gradu_gradv_tp(space, space, msh, problem_data.c_diff);
    mass_mat = op_u_v_tp(space, space, msh, problem_data.c_cap);
    rhs = op_f_v_time_tp(space, msh, problem_data, itime);
    
    % Apply Neumann boundary conditions
    for iside = nmnn_sides
        if (msh.ndim > 1)
            % Restrict the function handle to the specified side, in any dimension, gside = @(x,y) g(x,y,iside)
            gside = @(varargin) g(varargin{:},iside);
            dofs = space.boundary(iside).dofs;
            rhs(dofs) = rhs(dofs) + op_f_v_tp (space.boundary(iside), msh.boundary(iside), gside);
        else
            if (iside == 1)
                x = msh.breaks{1}(1);
            else
                x = msh.breaks{1}(end);
            end
            sp_side = space.boundary(iside);
            rhs(sp_side.dofs) = rhs(sp_side.dofs) + g(x,iside);
        end
    end
    
    % Apply Dirichlet boundary conditions
    if ~isempty(problem_data.h)
        u = zeros (space.ndof, 1);
        [u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (space, msh, h, drchlt_sides);
        u(drchlt_dofs) = u_drchlt;
        
        int_dofs = setdiff (1:space.ndof, drchlt_dofs);
        rhs(int_dofs) = rhs(int_dofs) - stiff_mat(int_dofs, drchlt_dofs)*u(drchlt_dofs) + mass_mat(int_dofs, int_dofs)*u_prev(int_dofs);
        
        % Solve the linear system
        delta_t = problem_data.time_discretization(itime)-problem_data.time_discretization(itime+1);
        lhs = mass_mat(int_dofs, int_dofs) - stiff_mat(int_dofs, int_dofs) * delta_t;
        u(int_dofs) =  lhs\ rhs(int_dofs);
    else
        rhs = rhs + mass_mat * u_prev;
        
        % Solve the linear system
        delta_t = problem_data.time_discretization(itime+1) - problem_data.time_discretization(itime+1);
        lhs = mass_mat - stiff_mat * delta_t;
        u =  lhs\ rhs;
        u_prev = u;
    end
    
end % END TIME INTEGRATION LOOP

end