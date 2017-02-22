% HSPACE_COARSEN: coarsen the hierarchical space, updating the fields of the object.
%
%   [hspace, Ccoarse] = hspace_coarsen (hspace, hmsh, funs_to_reactivate, removed_cells, reactivated_cell)
%
% INPUT:
%
%   hspace:             object representing the fine hierarchical space (see hierarchical_space)
%   hmsh:               object representing the hierarchical mesh, already coarsened (see hierarchical_mesh)
%   funs_to_reactivate: cell array with the indices, in the tensor product setting, of the functions to reactivate for each level
%   removed_cells:      cell array with the elements removed during coarsening, for each level
%   reactivated_cell:   cell array with the elements to be reactivated
%
% OUTPUT:
%
%   hspace:      object representing the coarsened hierarchical space (see hierarchical_space)
%   u_coarse:    coarsened dofs
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

function [hspace, u] = hspace_coarsen (hspace, hmsh, FtR, removed_cells, reactivated_cell)

boundary = ~isempty (hspace.boundary);

% Update active functions
if(nargout < 2)
    hspace = update_active_functions (hspace, hmsh, FtR, reactivated_cell);
else
    [hspace, u_coarse] = update_active_functions (hspace, hmsh, FtR, reactivated_cell, removed_cells);
end

% Update the matrices for changing basis
hspace.Csub = hspace_subdivision_matrix (hspace, hmsh);
u = zeros(hspace.ndof, 1);

for lev = hspace.nlevels:-1:1
    if ~isempty (u_coarse{lev})
        u = u + hspace.Csub{lev}'*u_coarse{lev};
    end
end


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
% function [hspace, u_coarse] = update_active_functions (hspace, hmsh, funs_to_reactivate, reactivated_cell)
%
% This function updates the active (hspace.active) and deactivated (hspace.deactivated) degrees of freedom,
% reactivating the functions in marked_funs.
% The function also updates hspace.nlevels, hspace.ndof and hspace.ndof_per_level
%
% Input:    hspace:                     the fine space, an object of the class hierarchical_space
%           hmsh:                       an object of the class hierarchical_mesh, already coarsened
%           funs_to_reactivate:         indices of active functions of level lev to be reactivated
%           reactivated_cell:           indices of cells of level lev to be reactivated
%
% Output:   hspace:    the coarsened space, an object of the class hierarchical_space
%           u_coarse:  coarsened_dofs
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

function [hspace, u_coarse] = update_active_functions (hspace, hmsh, funs_to_reactivate, reactivated_cell, removed_cells)

active = hspace.active;
deactivated = hspace.deactivated;
active_funs_supported = cell(1, hspace.nlevels); % cell array with the active functions of level l supported by the reactivated_cell
Bzr_ext_container = cell(1, hspace.nlevels);     % cell array to cache Bézier extraction matrices of each level
removed_funs = cell(1, hspace.nlevels);          % cell array eith functions to remove at each level

% initialize the cell array to store the coarsened dofs u_coarse
u_coarse = cell(1, hspace.nlevels);
ndof_per_level_prev = cellfun (@numel, hspace.active);
ndof_until_lev = sum (ndof_per_level_prev(1:hspace.nlevels-1));
ndof_above_lev = sum (ndof_per_level_prev(1:hspace.nlevels));
u_coarse{hspace.nlevels} = hspace.Csub{hspace.nlevels}(:,ndof_until_lev+1:ndof_above_lev)*hspace.dofs(ndof_until_lev+1:ndof_above_lev);

for lev = hspace.nlevels-1:-1:1
    if  ~isempty (reactivated_cell{lev})
        ndof_until_lev = sum (ndof_per_level_prev(1:lev-1));
        ndof_above_lev = sum (ndof_per_level_prev(1:lev));
        
        u_coarse{lev} = hspace.Csub{lev}(:,ndof_until_lev+1:ndof_above_lev)*hspace.dofs(ndof_until_lev+1:ndof_above_lev);
        active_funs_supported{lev} = intersect(sp_get_basis_functions (hspace.space_of_level(lev), hmsh.mesh_of_level(lev), reactivated_cell{lev}), active{lev});
        active{lev} = union (active{lev}, funs_to_reactivate{lev});
        deactivated{lev} = setdiff (deactivated{lev}, funs_to_reactivate{lev});
        
        removed_funs{lev+1} = intersect(sp_get_basis_functions(hspace.space_of_level(lev+1), hmsh.mesh_of_level(lev+1), removed_cells{lev+1}), hspace.active{lev+1});
        active{lev+1} = setdiff (active{lev+1}, removed_funs{lev+1});
        deactivated{lev+1} = union (deactivated{lev+1}, removed_funs{lev+1});
        u_coarse{lev+1}(removed_funs{lev+1}) = zeros(size(removed_funs{lev+1})); 
        
    end
end

if (~hspace.truncated)
    disp('ERROR: THB-spline are required for coarsening !!!');
end

% Computation of the matrix to pass from the original to the coarsened space
if (nargout == 2)
    % Bèzier extraction of the finest level
    knots = hspace.space_of_level(hspace.nlevels).knots;
    degree = hspace.space_of_level(hspace.nlevels).degree;
    Bzr_ext_container{hspace.nlevels} = bzrextr(knots{1}, degree(1));
    
    % loop over levels from n-1 to 0
    for lev = hspace.nlevels-1:-1:1
        if  ~isempty (reactivated_cell{lev})
            
            % Bèzier extraction of level lev
            knots = hspace.space_of_level(lev).knots;
            degree = hspace.space_of_level(lev).degree;
            Bzr_ext_container{lev} = bzrextr(knots{1}, degree(1));
            
            % express active_funs_supported as linear combination of funs of
            % level lev+1
            ndof_until_lev = sum (ndof_per_level_prev(1:lev));
            ndof_above_lev = sum (ndof_per_level_prev(1:lev+1));
            
            active_dofs_lc = hspace.Csub{lev+1}(:,active_funs_supported{lev})*hspace.dofs(active_funs_supported{lev});
            active_dofs_finelev = hspace.Csub{lev+1}(:,ndof_until_lev+1:ndof_above_lev)*hspace.dofs(ndof_until_lev+1:ndof_above_lev);
            active_funs_to_project = find(cat(1, active_dofs_lc, active_dofs_finelev));                                          % index of the funs to be projected
            cells_activated_in_proj = sp_get_cells(hspace.space_of_level(lev), hmsh.mesh_of_level(lev), active_funs_to_project); % cells to be projected
            
            u_coarse_temp = cell(1, numel(cells_activated_in_proj));
            
            hsource_el = 1/hmsh.mesh_of_level(lev+1).nel;
            htarget_el = 1/hmsh.mesh_of_level(lev).nel;
            
            % loop over elements to be reactivated
            for el = 1:numel(cells_activated_in_proj)
                % initialize u_coarse_temp of the element el
                u_coarse_temp{el} = zeros(size(active{lev}));
                % get indeces of functions to be activated in projection
                funs_projected = sp_get_basis_functions (hspace.space_of_level(lev), hmsh.mesh_of_level(lev), cells_activated_in_proj(el));
                % get children of the reactivated cell
                children_cells = hmsh_get_children(hmsh, lev, cells_activated_in_proj(el));
                % evaluate Bézier projector
                B_el_proj = bzrproj_el_hcoarse( hspace.space_of_level(lev).degree(1),...
                    Bzr_ext_container{lev+1}, inv(Bzr_ext_container{lev}(:,:,cells_activated_in_proj(el))), hsource_el, htarget_el, children_cells);
                % loop over children of target elements
                for child = 1:numel(children_cells)
                    % get indeces of functions with support on children cells to be
                    % projected
                    funs_supported = sp_get_basis_functions (hspace.space_of_level(lev+1), hmsh.mesh_of_level(lev+1), children_cells(child));
                    %project
                    u_coarse_temp{el}(funs_projected) = u_coarse_temp{el}(funs_projected) + B_el_proj(:,:,child)*(active_dofs_lc(funs_supported)+active_dofs_finelev(funs_supported));
                end
            end
            % end elements loop
            
            % smoothing of active functions with support on reactiveted cells
            funs_to_smooth = union(active_funs_supported{lev}, funs_to_reactivate{lev});
            
            % loop over active functions
            [~, funs_support_index] = sp_get_cells (hspace.space_of_level(lev), hmsh.mesh_of_level(lev), funs_to_smooth);

            smooth_dofs = zeros(size(funs_to_smooth));
            
            for i=1:numel(funs_to_smooth)
                for j=1:numel(funs_support_index{i})
                   funs_supported_smoothing = sp_get_basis_functions (hspace.space_of_level(lev+1), hmsh.mesh_of_level(lev+1), funs_support_index{i}(j));
                    for k=1:numel(funs_supported_smoothing)
                        if ~isempty (intersect(funs_supported_smoothing(k), funs_to_smooth(i)))
                            smooth_dofs(i) = smooth_dofs(i) + u_coarse_temp{funs_support_index{i}(j)}(funs_to_smooth(i))/numel(funs_support_index{i});
                        end
                    end
                end
            end
            % end smoothing
            
            % instert the new dofs into the hierarchical structure
            u_coarse{lev}(funs_to_smooth) = smooth_dofs;
            
        end
        
    end
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
    hspace.coeff_pou = u_coarse_temp * hspace.coeff_pou;
end

end
