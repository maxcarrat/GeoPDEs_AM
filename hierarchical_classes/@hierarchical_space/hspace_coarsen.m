% HSPACE_COARSEN: coarsen the hierarchical space, updating the fields of the object.
%
%   [hspace, u] = hspace_coarsen (hspace, hmsh, u, funs_to_reactivate, removed_cells)
%
% INPUT:
%
%   hspace:         object representing the fine hierarchical space (see hierarchical_space)
%   hmsh:           object representing the hierarchical mesh, already coarsened (see hierarchical_mesh)
%   u:              degrees of freedom to project onto the coarser mesh
%   funs_to_reactivate: cell array with the indices, in the tensor product setting, of the functions to reactivate for each level
%   removed_cells:  cell array with the elements removed during coarsening, for each level
%
% OUTPUT:
%
%   hspace:    object representing the coarsened hierarchical space (see hierarchical_space)
%   Ccoarse:   
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

function [hspace, u] = hspace_coarsen (hspace, hmsh, u, FtR, removed_cells)

boundary = ~isempty (hspace.boundary);

% Update of active functions
if (nargout == 2)
  [hspace, u] = update_active_functions (hspace, hmsh, u, FtR, removed_cells);
else
  hspace = update_active_functions (hspace, hmsh, u, FtR, removed_cells);
end

% % Update the matrices for changing basis
% hspace.Csub = hspace_subdivision_matrix (hspace, hmsh);
% 
% % Fill the information for the boundaries
% if (boundary)% && hmsh.ndim > 1)
%     Nf = cumsum ([0, hspace.ndof_per_level]);
%     for iside = 1:2*hmsh.ndim
%         if (hmsh.ndim > 1)
%             FtR_boundary = cell (size (FtR));
%             for lev = 1:numel (FtR)
%                 [~,~,FtR_boundary{lev}] = intersect (FtR{lev}, hspace.space_of_level(lev).boundary(iside).dofs);
%             end
%             cells_boundary = cell (size (removed_cells));
%             for lev = 1:numel (removed_cells)
%                 cells_boundary{lev} = get_boundary_indices (iside, hmsh.mesh_of_level(lev).nel_dir, removed_cells{lev});
%             end
%             hspace.boundary(iside) = hspace_coarsen (hspace.boundary(iside), hmsh.boundary(iside), FtR_boundary, cells_boundary);
%             
%             nlevels_aux = hspace.boundary(iside).nlevels;
%         elseif (hmsh.ndim == 1)
%             nlevels_aux = hspace.nlevels;
%         end
%         
%         dofs = [];
%         for lev = 1:nlevels_aux;
%             [~,iact] = intersect (hspace.active{lev}, hspace.space_of_level(lev).boundary(iside).dofs);
%             dofs = union (dofs, Nf(lev) + iact);
%         end
%         hspace.boundary(iside).dofs = dofs;
%     end
%     
% else
%     hspace.boundary = [];
% end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [hspace, u] = update_active_functions (hspace, hmsh, u, marked_funs, removed_cells)
%
% This function updates the active (hspace.active) and deactivated (hspace.deactivated) degrees of freedom,
% reactivating the functions in marked_funs.
% The function also updates hspace.nlevels, hspace.ndof and hspace.ndof_per_level
%
% Input:    hspace:    the fine space, an object of the class hierarchical_space
%           hmsh:      an object of the class hierarchical_mesh, already coarsened
%           u:         degrees of freedom to project onto the coarser mesh
%           marked_funs{lev}: indices of active functions of level lev to be reactivated
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

function [hspace, u] = update_active_functions (hspace, hmsh, u, funs_to_reactivate, removed_cells)

active = hspace.active;
ndof_per_level_cached = cellfun (@numel, active);
deactivated = hspace.deactivated;
removed_funs = cell(1, hspace.nlevels);
cell_to_reactivate = cell(1, hspace.nlevels);

for lev = hspace.nlevels:-1:2
    if (strcmpi (hspace.type, 'standard') && ~isempty (removed_cells{lev}))
        removed_funs{lev} = sp_get_basis_functions (hspace.space_of_level(lev), hmsh.mesh_of_level(lev), removed_cells{lev});
        cell_to_reactivate{lev-1} = sp_get_cells (hspace.space_of_level(lev-1), hmsh.mesh_of_level(lev-1), funs_to_reactivate{lev-1});
        active{lev} = setdiff (active{lev}, removed_funs{lev});
        active{lev-1} = union (active{lev-1}, funs_to_reactivate{lev-1});
        deactivated{lev-1} = setdiff (deactivated{lev-1}, funs_to_reactivate{lev-1});
        
    elseif (strcmpi (hspace.type, 'simplified') && ~isempty (funs_to_reactivate{lev-1}))
        active{lev-1} = union (active{lev-1}, funs_to_reactivate{lev-1});
        deactivated{lev-1} = setdiff (deactivated{lev-1}, funs_to_reactivate{lev-1});
        children = hspace_get_children (hspace, lev-1, funs_to_reactivate{lev-1});
        cell_to_reactivate{lev-1} = sp_get_cells (hspace.space_of_level(lev-1), hmsh.mesh_of_level(lev-1), funs_to_reactivate{lev-1});
        
        neighbors = sp_get_neighbors (hspace.space_of_level(lev-1), hmsh.mesh_of_level(lev-1), funs_to_reactivate{lev-1});
        deact_neighs = intersect (deactivated{lev-1}, neighbors);
        children = setdiff (children, hspace_get_children (hspace, lev-1, deact_neighs));
        active{lev} = setdiff (active{lev}, children);
    end
end

%% Bézier projection 
% project the fine level dofs to be coarsened onto the underlying coarser
% level

if (nargout == 2 || ~hspace.truncated)
    
    %HB-splines case
    
    ndof_per_level = cellfun (@numel, active);
    ndlev = hspace.ndof_per_level(hspace.nlevels);
    active_and_deact = union (active{hspace.nlevels}, deactivated{hspace.nlevels});
    [~,indices] = intersect (active_and_deact, hspace.active{hspace.nlevels});
    Id = sparse (numel(active_and_deact), ndlev);
    Id(indices,:) = speye (ndlev, ndlev);
    Ccoarse = Id;
    
    % update the hierarchical space
    hspace.active = active(1:hspace.nlevels);
    hspace.deactivated = deactivated(1:hspace.nlevels);
    hspace.ndof_per_level = cellfun (@numel, hspace.active);
    hspace.ndof = sum (hspace.ndof_per_level);
    
    % construct the sub-division operator for the coarse mesh
    Csub_coarsened = hspace_subdivision_matrix (hspace, hmsh, 'reduced');
    
    u_coarse = cell(1, hspace.nlevels);
    for lev = hspace.nlevels:-1:2
        % coarse only elements with marked cells
        if (~isempty(removed_cells{lev}))
            % get all the dofs of level lev
            u_ref = hspace.Csub{lev}*u(1:sum(ndof_per_level_cached(1:lev)));
            degree = hspace.space_of_level(hmsh.nlevels-1).degree;
            removed_funs_level = removed_funs{lev};
            cells_to_reactivate = cell_to_reactivate{lev-1};
            % obtain the dofs of the mesh @lev-1 projecting all the dofs of
            % lev onto it
            u_coarse{lev} = sp_coarse(hspace.space_of_level(lev), hspace.space_of_level(lev-1),...
                hmsh.mesh_of_level(lev), hmsh.mesh_of_level(lev-1), u_ref, removed_funs_level,...
                cells_to_reactivate, hmsh.nsub, degree, degree-1);
            % reconstruct the hierarchical space
            for lev_rec = lev:1:hspace.nlevels-1
                u = Csub_coarsened{lev_rec-1}'*u_coarse{lev_rec};
            end
%             [~,deact_indices] = intersect (active_and_deact, deactivated{lev});
%             [~,act_indices] = intersect (active_and_deact, active{lev});
%             active_and_deact = union (active{lev+1}, deactivated{lev+1});
%             
%             ndof_prev_levs = sum (ndof_per_level(1:lev-1));
%             ndof_until_lev = sum (ndof_per_level(1:lev));
%             
%             aux = sparse (ndof_until_lev + numel(active_and_deact), size(Ccoarse,2));
%             aux(1:ndof_prev_levs,:) = Ccoarse(1:ndof_prev_levs,:);
%             aux(ndof_prev_levs+(1:numel(active{lev})),:) = Ccoarse(ndof_prev_levs+act_indices,:);
%             aux(ndof_until_lev+(1:numel(active_and_deact)),:) = ...
%                 Cmat(active_and_deact,deactivated{lev}) * Ccoarse(ndof_prev_levs+deact_indices,:);
%             
%             ndlev = hspace.ndof_per_level(lev+1);
%             [~,indices] = intersect (active_and_deact, hspace.active{lev+1});
%             Id = sparse (numel(active_and_deact), ndlev);
%             Id(indices,:) = speye (ndlev, ndlev);
%             
%             Ccoarse = [aux, [sparse(ndof_until_lev,ndlev); Id]];
%             clear aux;
        end
    end
end

hspace.dofs = u;
hspace.Csub = Csub_coarsened;




end
