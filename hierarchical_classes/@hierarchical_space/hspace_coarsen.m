% HSPACE_COARSEN: coarsen the hierarchical space, updating the fields of the object.
%
%   [hspace, Ccoarse] = hspace_coarsen (hspace, hmsh, funs_to_reactivate, removed_cells)
%
% INPUT:
%
%   hspace:         object representing the fine hierarchical space (see hierarchical_space)
%   hmsh:           object representing the hierarchical mesh, already coarsened (see hierarchical_mesh)
%   funs_to_reactivate: cell array with the indices, in the tensor product setting, of the functions to reactivate for each level
%   removed_cells:  cell array with the elements removed during coarsening, for each level
%
% OUTPUT:
%
%   hspace:    object representing the coarsened hierarchical space (see hierarchical_space)
%   Ccoarse:   matrix to pass from the original to the coarsened space
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

function [hspace, Ccoarse] = hspace_coarsen (hspace, hmsh, FtR, removed_cells)

boundary = ~isempty (hspace.boundary);

% Update active functions
if(nargout < 2)
    hspace = update_active_functions (hspace, hmsh, FtR, removed_cells);
else
   [hspace, Ccoarse] = update_active_functions (hspace, hmsh, FtR, removed_cells);
end

% Update the matrices for changing basis
%% TODO make this function working for coarsening
hspace.Csub = hspace_subdivision_matrix (hspace, hmsh);

% Fill the information for the boundaries
if (boundary)% && hmsh.ndim > 1)
  Nf = cumsum ([0, hspace.ndof_per_level]);
  for iside = 1:2*hmsh.ndim
    if (hmsh.ndim > 1)
      FtR_boundary = cell (size (FtR));
      for lev = 1:numel (FtR)
        [~,~,FtR_boundary{lev}] = intersect (FtR{lev}, hspace.space_of_level(lev).boundary(iside).dofs);
      end
      cells_boundary = cell (size (removed_cells));
      for lev = 1:numel (removed_cells)
        cells_boundary{lev} = get_boundary_indices (iside, hmsh.mesh_of_level(lev).nel_dir, removed_cells{lev});
      end
      hspace.boundary(iside) = hspace_coarsen (hspace.boundary(iside), hmsh.boundary(iside), FtR_boundary, cells_boundary);

      nlevels_aux = hspace.boundary(iside).nlevels;
    elseif (hmsh.ndim == 1)
      nlevels_aux = hspace.nlevels;
    end
      
    dofs = [];
    for lev = 1:nlevels_aux;
      [~,iact] = intersect (hspace.active{lev}, hspace.space_of_level(lev).boundary(iside).dofs);
      dofs = union (dofs, Nf(lev) + iact);
    end
    hspace.boundary(iside).dofs = dofs;
  end
  
else
  hspace.boundary = [];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function hspace = update_active_functions (hspace, hmsh, funs_to_reactivate, removed_cells)
%
% This function updates the active (hspace.active) and deactivated (hspace.deactivated) degrees of freedom,
% reactivating the functions in marked_funs. 
% The function also updates hspace.nlevels, hspace.ndof and hspace.ndof_per_level
%
% Input:    hspace:    the fine space, an object of the class hierarchical_space
%           hmsh:      an object of the class hierarchical_mesh, already coarsened
%           funs_to_reactivate{lev}: indices of active functions of level lev to be reactivated
%           removed_cells{lev}: indices of cells of level lev that have been removed 
%
% Output:   hspace:    the coarsened space, an object of the class hierarchical_space
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

function [hspace, Ccoarse] = update_active_functions (hspace, hmsh, funs_to_reactivate, removed_cells)

active = hspace.active;
deactivated = hspace.deactivated;
removed_funs = cell(1, hspace.nlevels);

for lev = hspace.nlevels:-1:2
  if (strcmpi (hspace.type, 'standard') && ~isempty (removed_cells{lev}))
    removed_funs{lev} = sp_get_basis_functions (hspace.space_of_level(lev), hmsh.mesh_of_level(lev), removed_cells{lev});
    deactivated{lev} = intersect (active{lev}, removed_funs{lev});
    active{lev} = setdiff (active{lev}, removed_funs{lev});
    active{lev-1} = union (active{lev-1}, funs_to_reactivate{lev-1});
%     deactivated{lev-1} = setdiff (deactivated{lev-1}, funs_to_reactivate{lev-1});
    
  elseif (strcmpi (hspace.type, 'simplified') && ~isempty (funs_to_reactivate{lev-1}))
    active{lev-1} = union (active{lev-1}, funs_to_reactivate{lev-1});
    deactivated{lev-1} = setdiff (deactivated{lev-1}, funs_to_reactivate{lev-1});
    children = hspace_get_children (hspace, lev-1, funs_to_reactivate{lev-1});

    neighbors = sp_get_neighbors (hspace.space_of_level(lev-1), hmsh.mesh_of_level(lev-1), funs_to_reactivate{lev-1});
    deact_neighs = intersect (deactivated{lev-1}, neighbors);
    children = setdiff (children, hspace_get_children (hspace, lev-1, deact_neighs));
    active{lev} = setdiff (active{lev}, children);
  end
end

% Computation of the matrix to pass from the original to the coarsened space
if (nargout == 2 || ~hspace.truncated)
    
    %HB-splines case
    
    ndof_per_level = cellfun (@numel, active);
    ndof = sum(ndof_per_level(:));
    ndof_until_level = sum(ndof_per_level(1:end-1));
    ndlev = hspace.ndof_per_level(hspace.nlevels);
    active_and_deact = union (active{hspace.nlevels}, deactivated{hspace.nlevels});
    [~,indices] = intersect (hspace.active{hspace.nlevels}, active_and_deact);
    Id = sparse (ndof, ndlev);
    Id(ndof_until_level+indices,:) = speye (ndlev, ndlev);
    Ccoarse = Id;
    aux = sparse (ndof, size(Ccoarse,2));

    
    for lev = hspace.nlevels-1:-1:1
        % project the source space of level lev+1 onto the target space
        % of level lev
        if (~isempty(removed_cells{lev+1}))
            degree = hspace.space_of_level(lev).degree;
            [~, Proj] = sp_coarse (hspace.space_of_level(lev+1), hspace.space_of_level(lev), degree, degree-1);
            for idim=1:hmsh.ndim
                hspace.Proj{lev,idim} = Proj{idim};
            end
        end
%         Cmat = matrix_basis_change_coarsening (hspace, lev, deactivated{lev+1});
        Cmat = matrix_basis_change__ (hspace, lev+1);

        [~,removed_func_indeces] = intersect (active_and_deact, removed_funs{lev+1});
        
        active_and_deact = union (active{lev}, deactivated{lev});
%         [~,deact_indices] = intersect (active_and_deact, deactivated{lev});
        [~,act_indices] = intersect (active{lev}, hspace.active{lev});
        
        ndof_prev_levs = sum (ndof_per_level(1:lev));
        ndof_until_lev = sum (ndof_per_level(1:lev-1));
        
        ndlev = hspace.ndof_per_level(lev);

        aux(ndof_until_lev+act_indices,:) = Ccoarse(ndof_until_lev+act_indices,:);
%         aux(ndof_until_lev+(1:numel(active{lev})),:) = Ccoarse(ndof_until_lev+act_indices,:);
        aux(ndof_until_lev+(1:numel(funs_to_reactivate{lev})),:) = ...
            Cmat(funs_to_reactivate{lev},ndof_until_lev+removed_func_indeces) * Ccoarse(ndof_until_lev+removed_func_indeces,:);
        
        [~,indices] = intersect (active{lev}, hspace.active{lev});
        Id = sparse (ndof, ndlev);
        Id(ndof_until_lev+indices,:) = speye (numel(indices), ndlev);
        
%         ndlev_prev = sum(hspace.ndof_per_level(lev+1:end));
        Ccoarse = [Id aux];
        aux = Ccoarse;

    end
    
else
    error('THBsplines not implemented yet');
end

hspace.active = active(1:hspace.nlevels);
hspace.deactivated = deactivated(1:hspace.nlevels);
hspace.ndof_per_level = cellfun (@numel, hspace.active);
hspace.ndof = sum (hspace.ndof_per_level);

hspace.active = active(1:hspace.nlevels);
hspace.deactivated = deactivated(1:hspace.nlevels);
hspace.ndof_per_level = cellfun (@numel, hspace.active);
hspace.ndof = sum (hspace.ndof_per_level);

if (hspace.truncated)
  hspace.coeff_pou = ones (hspace.ndof, 1);
else
  hspace.coeff_pou = Ccoarse * hspace.coeff_pou;
end

end
