function [ marked_element_index, num_marked_element ] = refine_toward_source( hmsh, hspace, time_step, adaptivity_data, problem_data )
%REFINE_TOWARD_SOURCE Generate a refined mesh around a given source defined
%vy a vertex and a radius.

msh_finer_lev = hmsh.mesh_of_level(end);
vertex = problem_data.path(time_step,:);

% MARK ELEMENT OF THE FINER LEVEL CUTTED OR INCLUDED BY THE SOURCE REGION
el_dir = zeros(3, 1);
marked_element_index = cell(hmsh.nlevels, 1);
for lev=1:hmsh.nlevels
    % get local coordinates
    local_vertex_left = mapGlobalToLocal(vertex - adaptivity_data.radius, msh_finer_lev);
    local_vertex_right = mapGlobalToLocal(vertex + adaptivity_data.radius, msh_finer_lev);
    % LOOP OVER ACTIVE DOFS OF THE LEVEL
    for el = 1:hmsh.nel_per_level(lev)
        isMarked=0;
        [el_dir(1), el_dir(2), el_dir(3)] = ind2sub(hmsh.mesh_of_level(lev).nel_dir, hmsh.active{lev}(el));
        % LOOP OVER DIRECTIONS
        for idir = 1:hmsh.rdim
            degree_dir = hspace.space_of_level(lev).degree(idir);
            knots_dir = hspace.space_of_level(lev).knots{idir};
            % check if vertex is inside the knot-span in ith direction
            if (((local_vertex_left(idir) > knots_dir(el_dir(idir)+degree_dir) &&... the source is contained within the knotspan
                    local_vertex_right(idir) <= knots_dir(el_dir(idir)+degree_dir+1)) ||...
                    (local_vertex_left(idir) < knots_dir(el_dir(idir)+degree_dir) &&... the knotspan is within the source
                    local_vertex_right(idir) >= knots_dir(el_dir(idir)+degree_dir+1)) ||...
                    (local_vertex_left(idir) < knots_dir(el_dir(idir)+degree_dir) &&... the knotspan is cutted by the source
                    local_vertex_right(idir) >= knots_dir(el_dir(idir)+degree_dir)) ||...
                    (local_vertex_left(idir) < knots_dir(el_dir(idir)+degree_dir+1) &&... the knotspan is cutted by the source
                    local_vertex_right(idir) >= knots_dir(el_dir(idir)+degree_dir+1))) && ...  is not the max level
                    lev <= adaptivity_data.max_level )
                isMarked = 1;
            else
                isMarked = 0;
                break;
            end
        end % END DIRECTION LOOP
        if isMarked
            % add element to the marked list
            marked_element_index{lev} = [marked_element_index{lev} hmsh.active{lev}(el)];
            if adaptivity_data.mark_neighbours
                % check if its neighbours are active
                global_element_index = sub2ind(hmsh.mesh_of_level(lev).nel_dir, el_dir(1), el_dir(2), el_dir(3));
                [neighbours, flag]  = hmsh_get_neighbours( hmsh, lev, global_element_index );
                n_neighbours = numel(neighbours);
                % if neighbours of the same level and the level below are
                % active, add the element in the marked list.
                % loop over neighbours
                for i=1:n_neighbours
                    % check if neighbour is of lower level...
                    if (flag(i)>1)
                        % ... and add it to the marked list
                        marked_element_index{lev} = [marked_element_index{lev} neighbours(i)];
                    end
                end % end loop over neighbours
            end
        end
    end % END ELEMENTS LOOP
end % END LEVELS LOOP
num_marked_element = sum(cellfun(@numel, marked_element_index));

end


function local_coords = mapGlobalToLocal(global_coords, msh)
local_origin = zeros(size(global_coords));
local_length = ones(size(global_coords));

local_coords =  (global_coords-msh.map(local_origin)')./(msh.map(local_length)'-msh.map(local_origin)');
local_coords = local_coords';
end

