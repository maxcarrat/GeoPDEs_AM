%HMSH_GET_NEIGHBOURS Get the elements in the neighbourhood of a given
%element. The neighbourhood of an element includes all the elements of the
%lower and the same level, which partially share either an edge (2D case) or a face
%(3D case) with the element itself.
%
%     [neighbours, flag] = hmsh_get_neighbours (hmsh, lev, ind)
%
%
% INPUT:
%
%     hmsh: the hierarchical mesh (see hierarchical_mesh)
%     lev:  level of the cells to subdivide
%     ind:  indices of the cells in the Cartesian grid
%
% OUTPUT:
%
%     neighbours: active neighbouring elements, numebered as in figure:
%
%                      n_top (6)
%                        +   n_rear (3)
%                        +  +
%     (1)  n_left   + + ind + + n_right (2)
%                      + +
%           (4) n_front  +
%                     n_bottom (5)
%
%     index (1 x neighbour cell array): if active and of the same level(2), if not
%                                       active(0), if active and of a lower
%                                       level(1), if active and of higher level(3).
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
function  [neighbours, index]  = hmsh_get_neighbours( hmsh, lev, ind )
neighbours = [];
index = [];
% get neighbours of the same level
neighbours_x_dir = [ind-1, ind+1];
neighbours_y_dir = [ind-hmsh.mesh_of_level(lev).nel_dir(1), ind+hmsh.mesh_of_level(lev).nel_dir(1)];
neighbours_z_dir = [ind-hmsh.mesh_of_level(lev).nel_dir(1)*ind-hmsh.mesh_of_level(lev).nel_dir(2), ind-hmsh.mesh_of_level(lev).nel_dir(1)*ind-hmsh.mesh_of_level(lev).nel_dir(2)];
% fill an auxiliary list of elements
aux = [neighbours_x_dir neighbours_y_dir neighbours_z_dir];
% loop over auxiliary indeces
for i=1:numel(aux)
    % check if neighbours are active
    if ~isempty(intersect(aux(i), hmsh.active{lev}))
        neighbours = [neighbours aux(i)];
        index = [index 2];
    else
        isMarked = 0;
        if (lev > 1 && aux(i) > 0 && aux(i)<hmsh.mesh_of_level(lev).nel)
            [parent, ~] = hmsh_get_parent(hmsh, lev, aux(i));
            % check if parents are active
            if ~isempty(intersect(parent, hmsh.active{lev-1}))
                neighbours = [neighbours aux(i)];
                index = [index 1];
                isMarked = 1;
            else
                neighbours = [neighbours aux(i)];
                index = [index 0];
            end
        elseif(lev < hmsh.nlevels && aux(i) > 0 && aux(i)<hmsh.mesh_of_level(lev).nel && ~isMarked)
            [children, ~] = hmsh_get_children(hmsh, lev, aux(i));
            if ~isempty(intersect(children, hmsh.active{lev+1}))
                neighbours = [neighbours aux(i)];
                index = [index 3];
            else
                neighbours = [neighbours aux(i)];
                index = [index 0];
            end
        end
    end % end if
end % end for loop

end

