% OP_F_V_TIME_TP: assemble the right-hand side time-dependent vector r(t) = [r(i)], with  r(i) = (f(t), v_i), exploiting the tensor product structure.
%
%   rhs = op_f_v_time_tp (spv, msh, coeff);
%
% INPUT:
%
%   spv:   object representing the function space (see sp_scalar)
%   msh:   object defining the domain partition and the quadrature rule (see msh_cartesian)
%   coeff: function handle to compute the source function
%
% OUTPUT:
%
%   rhs: assembled right-hand side
%
% Copyright (C) 2017 Massimo Carraturo
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

function rhs = op_f_v_time_tp (space, msh, problem_data, time_step)

rhs = zeros (space.ndof, 1);

for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col  = sp_evaluate_col (space, msh_col);
    x = cell (msh.rdim, 1);
    path = cell (msh.rdim, 1);
    for idim = 1:msh.rdim
        x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
        path{idim} = ones(size(x{idim})).*problem_data.path(time_step, idim);
    end
    
    rhs = rhs + op_f_v (sp_col, msh_col, problem_data.f(x{:}, path{:}));
end

end
