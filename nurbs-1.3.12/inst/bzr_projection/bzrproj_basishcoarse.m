function B = bzrproj_basishcoarse( p, C_source, C_target, taget_el_bound_vec )
%BZRPROJ_BASISHCOARSE: global Bézier projection operator (h-coarsening)
% This projection operator projects the basis functions of the source (fine) 
% mesh onto the ones of the target (coarse) mesh.
%
% Calling Sequence:
%
%   B = bzrproj_basishcoarse( p, C_source, C_target, taget_el_bound_vec )
%
%    INPUT:
%      p                        - polynomial degree
%      C_source                 - Bézier extraction operator of the source mesh
%      C_target                 - Bézier extraction operator of the source mesh
%      taget_el_bound           - target elements boundaries [nsub x 2] 
%
%    OUTPUT:
%
%      B - global Bézier projection operator [ndof_target x ndof_source]
%
%   The function is implemented for the univariate case and uniform knot spans.
%   Smoothing coeffients are taken as: w_A = 1/n_A_el, where n_A_el is the number
%   of supporting elements of the Ath function.
%
%   Adapted from "Bézier projection: A unified approach for local
%   projection and quadrature-free refinement and coarsening of NURBS and
%   T-splines with particular application to isogeometric design and
%   analysis", D.C. Thomas et al., Comput. Methods Appl. Mech. Engrg. 284
%   (2015)55-105
%
%   Bézier projection hcoarsening of basis functions. adapts the original
%   implementation given only for control values projection, using the
%   formulas given in:
%   "Isogeometric finite element data structures based on Bézier extraction
%   of NURBS", Border et al., Int. J. Num. Meth. Engrg. (2011) vol. 87 pp.
%   15-47.
%
% Copyright (C) 2017 Massimo Carraturo
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

% evaluate source and target parameters
ndofs_source = size(C_source, 3)+(p-1);         % number of functions of the source mesh
ndofs_target = size(C_target, 3)+(p-1);         % number of functions of the target mesh
nel_source = ndofs_source - p;                  % number of elements of the source mesh
nel_target = ndofs_target - p;                  % number of elements of the target mesh
h_source = 1/nel_source;                        % element size of the source mesh
h_target = 1/nel_target;                        % element size of the target mesh
n_sub = nel_source/nel_target;                  % number of source elements per target element

% initialize global Bézier projection operator
B = sparse(ndofs_target, ndofs_source);  

% index source elements vector
index_source_element_vec = linspace(1,ndofs_source,ndofs_source);

% build up matrix to convert Bézier coefficients of the source Bernstein
% basis into coefficients of the target Bernstein basis

% loop over target elements and assembly the corresponding Bézier
% projection operator
for i=1:nel_target
    % index source elements
    start_index = i + (i-1) * (n_sub-1) ;
    end_index = start_index + n_sub - 1;
    % projection operator of the ith target element
    B_target = bzrproj_el_basishcoarse( p, C_source, C_target(:,:,i), h_source, h_target,...
        taget_el_bound_vec(:,:), index_source_element_vec(start_index:end_index));
    % assembly the operators of each sub-element 
    B_target_el = zeros(p+1, n_sub + p);
    for j=1:n_sub
        B_target_el(:,j:j+p) = B_target_el(:,j:j+p) + B_target(:,:,j);
    end
    % assembly global projection operator 
    B(i:i+p, i + (n_sub-1)*(i-1) : (i-1) + p + n_sub + (n_sub-1)*(i-1)) =...
        B(i:i+p, i + (n_sub-1)*(i-1) : (i-1) + p + n_sub + (n_sub-1)*(i-1)) + B_target_el;
end

%% Smoothing ==============================================================
%
%     Weights for the Smoothing are tooken from "Convergence of an
%     efficient local least-squares fitting metho for bases with compact
%     support"; Govindjee, S., Strain, J., Mitchell T. J. and Taylor R. L.;
%     Comput. Methods. Appl. Mech. Engrg. 213-216 (2012) 84-92.
%
%

% evaluate weights for smoothing
weights = zeros(ndofs_target, 1);
for i=1:ndofs_target
    if (i < p + 1)
        weights(i) = 1/i;
    elseif (i > ndofs_target - p)
        weights(i) = 1/(ndofs_target-i+1);
    else
        weights(i) = 1/(p+1);
    end
    B(i,:) = weights(i)*B(i,:);
end

end

function [ B ] = bzrproj_el_basishcoarse( p, C_source, C_target, h_source, h_target, taget_el_bound_vec, index_source_element_vec )
%BZRPROJ_EL_HCOARSE:  Bézier projection operator (h-coarsening)
% This projection operator project n elements of a source mesh onto
% an element of a target mesh. The resulting control values have to be 
% smoothed by means of the weights.
%
% Calling Sequence:
%
%   B = bzrproj_hcoarse( p, C_source, R_target, h_source, h_target, taget_el_bound_vec, index_source_element_vec )
%
%    INPUT:
%      p                        - polynomial degree
%      C_source                 - Bézier extraction operator of the source mesh
%      C_target                 - Bézier extraction operator of the target mesh
%      h_source                 - element size of the source mesh
%      h_target                 - element size of the target mesh
%      taget_el_bound_vec       - vector of target element boundaries
%      index_source_element_vec - vector of index source elements
%
%    OUTPUT:
%
%      B - Bézier projection operator [p+1 x p+1 x n+1]
%
%   Adapted from "Bézier projection: A unified approach for local
%   projection and quadrature-free refinement and coarsening of NURBS and
%   T-splines with particular application to isogeometric design and
%   analysis", D.C. Thomas et al., Comput. Methods Appl. Mech. Engrg. 284
%   (2015)55-105
%
%   Algorithm 4.14 pp. 94-95
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

% initialize Bézier projection operator
n = numel(index_source_element_vec);
B = zeros(p+1, p+1, n);

% build up matrix to convert Bézier coefficients of the source Bernstein
% basis into coefficients of the target Bernstein basis

% Gramian matrix and inverse
G = grmn(p);
G_inv = grmninv(p);
% source elements intervals in parent coordinates
xi = linspace(-1, 1, n+1);
% elements volume rate
phi = h_source/h_target;

% loop over elements of the source mesh
for i=1:n
    % transformation matrix from [-1, 1] to taget element boundaries
%     A_inv = inv(brnsttransf(taget_el_bound_vec(i,1), taget_el_bound_vec(i,2), p));
	A = brnsttransf(xi(i), xi(i+1), p);
    % projection operator
    B(:,:,i) = C_target(:,:) * phi*G_inv*A'*G * inv(C_source(:,:,index_source_element_vec(i)));
end

end
