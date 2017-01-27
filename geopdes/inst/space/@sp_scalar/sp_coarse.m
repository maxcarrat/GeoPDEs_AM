% SP_COARSE: returns the projected degrees of freedom from a source (fine)
% space onto a target (coarse) space.
%
%     [Proj] = sp_refine sp_source, sp_target, msh_source, msh_target, u, removed_funs, cells_to_reactivate, nsub, degree, regularity);
%
% INPUTS:
%
%     sp_source:            the fine space, an object of the sp_scalar class (see sp_scalar)
%     sp_target:            the coarse space, an object of the sp_scalar class (see sp_scalar)
%     msh_source:           an object of the msh_cartesian class, the fine mesh (see msh_cartesian)
%     msh_target:           an object of the msh_cartesian class, the coarse mesh (see msh_cartesian)
%     u:                    vector of degrees of freedom to project
%     removed_funs:         list of dofs of the source mesh to remove
%     cells_to_reactivate:  list of elements to reactivate
%     nsub:                 number of uniform merging to apply on each knot span, and for each direction
%     degree:               degree of the fine space, and for each direction
%     regularity:           regularity for the new space, and for each direction
%
% OUTPUT:
%
%     Proj:    the coefficients relating 1D splines of the coarse and the fine spaces
%
%     Weights for the Smoothing are tooken from "Convergence of an
%     efficient local least-squares fitting metho for bases with compact
%     support"; Govindjee, S., Strain, J., Mitchell T. J. and Taylor R. L.;
%     Comput. Methods. Appl. Mech. Engrg. 213-216 (2012) 84-92.
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

function [sp_coarse, Proj] = sp_coarse (sp_source, sp_target, msh_source, msh_target, u, removed_funs, cells_to_reactivate, nsub, degree, regularity)

if (any (degree < sp_target.degree))
    error ('The given degree should be greater or equal than the degree of the coarse space')
end

if (any (regularity > degree - 1))
    error ('The regularity cannot be greater than degree-1')
end

Proj = cell (1, msh_target.ndim);
if (strcmpi (sp_target.space_type, 'spline'))
    if (all (sp_target.degree == degree))
        
        proj_temp_vec = cell(1, msh_target.nel);
        % initialize Bézier extraction operators of the source and target
        % mesh
        C_source = cell(1, msh_target.ndim);
        C_target = cell(1, msh_target.ndim);
        R_target = cell(1, msh_target.ndim);
        % element volume considering uniform knot vectors in each
        % direction
        h_coarse = 1;
        h_fine = 1;
        el_bounds = cell(msh_target.ndim);
        % loop over each dimension
        for idim = 1:msh_target.ndim
            h_coarse = h_coarse * 1/(numel(sp_target.knots{idim}-2*degree(idim)));             % coarse element size
            h_fine = h_fine * 1/(numel(sp_source.knots{idim}-2*degree(idim)));                 % fine element size
            % evaluate fine element parametric bounds
            nsub_dir = nsub(idim);
            el_bounds{idim} = eval_param_bounds(1/(numel(sp_target.knots{idim})-2*degree(idim)), 1/(numel(sp_source.knots{idim})-2*degree(idim)), nsub_dir);
            % source mesh Bézier estraction operator
            C_source{idim} = bzrextr(sp_source.knots{idim}, degree(idim));
            % target mesh Bézier estraction operator
            C_target{idim} = bzrextr(sp_target.knots{idim}, degree(idim));
        end
        
        % loop over coarse mesh elements in idim
        for icell_global = 1:msh_target.nel
             %% Bézier projection ------------------------------------------
            % project the dofs vector u onto the coarse mesh element
            icell_dir = cell(1, msh_target.ndim);                          % ith cell index in each direction
            sub_el_indeces = cell(1, msh_target.ndim);
            nel_dir = cell(1, msh_target.ndim);                            % number of target elements in each parametric direction
            nsubel_dir = cell(1, msh_source.ndim);                         % number of source elements in each parametric direction
            nsub_dir = nsub(1);
            nel_dir{1} = (sp_target.ndof_dir(1)-degree(1));
            nsubel_dir{1} = (sp_source.ndof_dir(1)-degree(1));
            icell_dir{1} = nel_dir{1} - (ceil(icell_global/nel_dir{1})*nel_dir{1}-icell_global);
            % evaluate subelements indeces
            start_bzr_index = 1 + nsub_dir*(icell_dir{1}-1);
            end_bzr_index = start_bzr_index + (nsub_dir-1);
            sub_el_indeces{1} =  linspace(start_bzr_index, end_bzr_index, nsub_dir);
            % Spline reconstruction operator of target mesh
            R_target{1} = inv(C_target{1}(:,:,icell_dir{1}));
            % count the elements in each directions
            nel_count = 1;
            
            % loop over each dimension
            for idim = 2:msh_target.ndim
                nsub_dir = nsub(idim);
                nel_dir{idim} = (sp_target.ndof_dir(idim)-degree(idim));
                nsubel_dir{idim} = (sp_source.ndof_dir(idim)-degree(idim));
                icell_dir{idim} = ceil(icell_global/(nel_dir{idim-1}*nel_count)-1.0e-18);
                % evaluate subelements indeces
                start_bzr_index = nsub_dir*(icell_dir{idim}-1) + 1;
                end_bzr_index = start_bzr_index + (nsub_dir-1);
                sub_el_indeces{idim} =  linspace(start_bzr_index, end_bzr_index, nsub_dir);
                % Spline reconstruction operator of target mesh
                R_target{idim} = inv(C_target{idim}(:,:,icell_dir{idim}));
                % update elements counter
                nel_count = nel_dir{idim-1};
            end % end loop over dimensions
            
            iter = 0;
            % loop over sub-elements:
            for idim = 1:msh_target.ndim
                for isub=1:nsub(idim)
                    
                    % 1) find the global index ogf the sub-cell to project
                    if msh_target.ndim == 1
                        isubcell_global = sub_el_indeces{idim}(isub);
                    elseif msh_target.ndim == 2
                        isubcell_global = sub_el_indeces{idim}(isub) + (sub_el_indeces{idim}(isub)-1)*nsubel_dir{idim};
                    else
                        isubcell_global = sub_el_indeces{idim}(isub) + (sub_el_indeces{idim}(isub)-1)*nsubel_dir{idim} + (sub_el_indeces{idim}(isub)-1)*nsubel_dir{idim}*nsubel_dir{idim};
                    end
                    
                    % 2) extract the function indeces of the ith sub-cell
                    isubfncts = sp_get_basis_functions (sp_source, msh_source, isubcell_global);
                end % end loop ovre sub-elements
                
            end % end loop over each parametric direction     
            
            % 3) project onto the target mesh
            % evaluate the tensor product Bézier operator
            B_tp = bzrproj_hcoarse_tp(degree, C_source, R_target, msh_target.ndim, isub, el_bounds, sub_el_indeces);
            % project the fine mesh dofs onto the coarse mesh
            %                         if (iter == 0)
            proj_temp_vec{icell_global} = ...
                h_fine/h_coarse * B_tp(:,:)*u(isubfncts);
            %                         else
            %                             proj_temp_vec{icell_global} = proj_temp_vec{icell_global} +...
            %                                 h_fine/h_coarse * B_tp(:,:)*u(isubfncts);
            %                         end
            
            iter = iter +1;
        end % end loop over cells to remove
        
        %% Smoothing ------------------------------------------------------
        % loop over degrees of freedom of target space 
        Proj = zeros(sp_target.ndof, 1);
        ndof_support = (degree(1)^msh_target.ndim);

        for idof=1:sp_target.ndof
            % evaluate function supports
            support_indeces = sp_get_cells (sp_target, msh_target, idof);
            % smoothing spline dofs
            for isupport=1:numel(support_indeces)
                if(~isempty(proj_temp_vec{isupport}))
                    for idof_support = 1:ndof_support
                        Proj(idof) = Proj(idof) + 1/numel(support_indeces) * proj_temp_vec{support_indeces(isupport)}(idof_support);
                    end
                end
                
            end % end loop over supports
            
        end % end loop over target dofs
        
    else
        error ('For the Proj matrix, the case of degree elevation is not implemented yet')
    end
elseif (strcmpi (space.space_type, 'nurbs'))
    % Overkilling to compute the new weights, but does not change that much
    coefs = zeros ([4, size(space.weights)]);
    coefs(4,:,:,:) = space.weights;
    
    if (numel (space.knots) ~= 1)
        aux_nurbs = nrbdegelev (nrbmak (coefs, space.knots), degree - space.degree);
    else
        aux_nurbs = nrbdegelev (nrbmak (coefs, space.knots{1}), degree - space.degree);
    end
    
    [knots,~,new_knots] = kntrefine (aux_nurbs.knots, nsub_dir-1, aux_nurbs.order-1, regularity);
    aux_nurbs = nrbkntins (aux_nurbs, new_knots);
    weights = reshape (aux_nurbs.coefs(4,:,:,:), [aux_nurbs.number, 1]); % The extra 1 makes things work in any dimension
    
    sp_fine = sp_nurbs (knots, degree, weights, msh_target);
    
    if (nargout == 2)
        if (all (sp_fine.degree == space.degree))
            for idim = 1:msh_target.ndim
                knt_coarse = space.knots{idim};
                knt_fine = sp_fine.knots{idim};
                Proj{idim} = basiskntins (degree(idim), knt_coarse, knt_fine);
            end
        else
            error ('For the Proj matrix, the case of degree elevation is not implemented yet')
        end
    end
end
end

function el_bounds = eval_param_bounds(h_coarse, h_fine, nsub)
% EVAL_PARM_BOUNDS: calculate the upper and lower bounds of each sub-element
% to be coarsened to the parametric coarse space [-1, 1], in case of uniform
% knot vector.
%
%     [el_bounds] = eval_param_bounds(h_coarse, h_fine, nsub)
%
% INPUTS:
%
%     h_coarse:      coarse element parametric size
%     h_fine:        fine element parametric size
%     nsub:          number of uniform merging to apply on each knot span, and for each direction
%
% OUTPUT:
%
%     el_bounds:    mapped bounds of size [nsub x 2]
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

el_bounds = zeros(nsub, 2);
for isub=1:nsub
    el_bounds(isub, 1) = -2*(isub-1)*h_fine/h_fine - 1;
    el_bounds(isub, 2) = 2*h_coarse/h_fine - 2*(isub-1)*h_fine/h_fine - 1;
end


%! test:
%! el_bounds = eval_param_bounds(1/2, 1/4, 2)
%! if(isequal(el_bounds,[[-1 3];[-3 1]]))
%!     disp('Success !!!');
%! else
%!     disp('Error !!!');
%! end

end


function support_index_vec = eval_support_index(p, ndofs, A)
% EVAL_SUPPORT_INDEX: evaluate the element support indeces of the Ath function.
%
%     [support_index_vec] = eval_support_index(p, ndof, index_dof)
%
%    INPUT:
%
%      p                    - polynomial degree
%      ndof                 - number of degrees of freedom
%      index_dof            - index of actual degree of freedom
%
% OUTPUT:
%
%     support_index_vec:    vector of support indeces of the Ath dofs
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

if (A < p + 1 )
    support_index_vec = linspace(1, A, A);
elseif (A > ndofs - (p+1))
    support_index_vec = linspace(A - p, ndofs - p, ndofs-A+1);
else
    support_index_vec = linspace(A - p, A, p+1);
end

end
