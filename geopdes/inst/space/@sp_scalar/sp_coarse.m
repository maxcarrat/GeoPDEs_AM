% SP_COARSE: returns the projected degrees of freedom from a source (fine)
% space onto a target (coarse) space.
%
%     [sp_coarse, Proj] = sp_coarse(sp_source, sp_target, nsub, degree, regularity);
%
% INPUTS:
%
%     sp_source:            the fine space, an object of the sp_scalar class (see sp_scalar)
%     sp_target:            the coarse space, an object of the sp_scalar class (see sp_scalar)
%     degree:               degree of the fine space, and for each direction
%     regularity:           regularity for the new space, and for each direction
%
% OUTPUT:
%
%     Proj:    the coefficients relating 1D splines of the coarse and the fine spaces
%
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

function [sp_coarse, Proj] = sp_coarse (sp_source, sp_target, degree, regularity)

if (any (degree < sp_target.degree))
    error ('The given degree should be greater or equal than the degree of the coarse space')
end

if (any (regularity > degree - 1))
    error ('The regularity cannot be greater than degree-1')
end

ndim = numel(degree);
Proj = cell (1, ndim);
if (strcmpi (sp_target.space_type, 'spline'))
    sp_coarse = sp_target;
    
    if (nargout == 2)
        if (all (sp_source.degree == sp_coarse.degree))
            for idim = 1:ndim
                % evaluate the Bézier extraction operators of the source and
                % target mesh
                C_target = bzrextr(sp_target.knots{idim}, degree(idim));
                C_source = bzrextr(sp_source.knots{idim}, degree(idim));
                % evaluate Bézier projection operator
                Proj{idim} = bzrproj_basishcoarse( degree(idim), C_source, C_target, cat(1, [-3 1], [-1 3]) );
            end
        else
            error ('For the Proj matrix, the case of degree elevation is not implemented yet')
        end
    end
    
elseif (strcmpi (sp_target.space_type, 'nurbs'))
    % Overkilling to compute the new weights, but does not change that much
    coefs = zeros ([4, size(sp_target.weights)]);
    coefs(4,:,:,:) = sp_target.weights;
    if (numel (sp_target.knots) ~= 1)
        aux_nurbs = nrbdegelev (nrbmak (coefs, sp_target.knots), degree - space.degree);
    else
        aux_nurbs = nrbdegelev (nrbmak (coefs, sp_target.knots{1}), degree - space.degree);
    end
    
    aux_nurbs = nrbkntins (aux_nurbs, new_knots);    
    sp_coarse = sp_target;
    
    if (nargout == 2)
        if (all (sp_coarse.degree == space_source.degree))
            for idim = 1:msh.ndim
                % evaluate the Bézier extraction operators of the source and
                % target mesh
                C_target = bzrextr(sp_target.knots{idim}, degree(idim));
                C_source = bzrextr(sp_source.knots{idim}, degree(idim));
                % evaluate Bézier projection operator
                Proj{idim} = bzrproj_basishcoarse( degree(idim), C_source, C_target, taget_el_bound_vec );
            end
        else
            error ('For the Proj matrix, the case of degree elevation is not implemented yet')
        end
    end
end

end