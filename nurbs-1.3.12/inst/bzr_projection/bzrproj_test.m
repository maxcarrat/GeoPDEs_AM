%% BZRPROJ_TEST
%   From "Bézier projection: A unified approach for local
%   projection and quadrature-free refinement and coarsening of NURBS and
%   T-splines with particular application to isogeometric design and
%   analysis", D.C. Thomas et al., Comput. Methods Appl. Mech. Engrg. 284
%   (2015)55-105
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

clc;
clear;
close all

%% Refine curve------------------------------------------------------------
% curve parameters
p = 2;
controlPoints = [[-0.1 -0.3 -0.15 0.7 0.9 1];[0.2 0.4 1.0 0.9 0.75 0.4]];
source_knots = [0 0 0 1/4 1/2 3/4 1 1 1];

% plot
crv = nrbmak(controlPoints, source_knots);
figure(1); 
nrbplot(crv, 200);
hold on;
grid on;

%% Bézier projection element level-----------------------------------------
% construct the Bézier extraction operator for the source and the target
% mesh

target_knots = [0 0 0 1/2 3/4 1 1 1];
refined_knots = [0 0 0 1/8 1/4 3/8 1/2 5/8 3/4 7/8 1 1 1];
C = bzrextr(source_knots, p);
C_coarse = bzrextr(target_knots, p);
C_refined = bzrextr(refined_knots, p);
complete_contolpoints = cat(2, [0 0 0 0 0]', [1 1 1 1 1]');

% project control values of the source mesh onto the target mesh
B_1 =  bzrproj_el_hcoarse( p, C, inv(C_coarse(:,:,1)), 1/4, 1/2, [1 2] );
B_2 =  bzrproj_el_hcoarse( p, C, inv(C_coarse(:,:,2)), 1/4, 1/4, 3 );
B_3 =  bzrproj_el_hcoarse( p, C, inv(C_coarse(:,:,3)), 1/4, 1/4, 4 );

% new Bézier elements control points
new_control_point_e1 = B_1(:,:,1) * crv.coefs(1:2, 1:3)' + B_1(:,:,2) * crv.coefs(1:2, 2:4)';
new_control_point_e2 = B_2(:,:,1) * crv.coefs(1:2, 3:5)';
new_control_point_e3 = B_3(:,:,1) * crv.coefs(1:2, 4:6)';

% smoothing using local least square weights
new_control_points(1,:) = new_control_point_e1(1,:);
new_control_points(2,:) = 1/2 * new_control_point_e1(2,:) + 1/2 * new_control_point_e2(1,:);
new_control_points(3,:) = 1/3 * new_control_point_e1(3,:) + 1/3 * new_control_point_e2(2,:) + 1/3 * new_control_point_e3(1,:);
new_control_points(4,:) = 1/2 * new_control_point_e2(3,:) + 1/2 * new_control_point_e3(2,:);
new_control_points(5,:) = new_control_point_e3(3,:);


% make B-spline
pcrv = nrbmak((cat(2, new_control_points, complete_contolpoints))', target_knots);

nrbplot(pcrv, 200);

%% Bézier projection global------------------------------------------------
% construct the Bézier extraction operator for the source and the target
% mesh

% h-coarsening
target_knots = [0 0 0 1/4 1/2 1 1 1];
C_coarse = bzrextr(target_knots, p);
complete_contolpoints = cat(2, [0 0 0 0 0]', [1 1 1 1 1]');

% project control values of the source mesh onto the target mesh
B_coarse =  bzrproj_hcoarse( p, C, C_coarse, 2);

new_control_points = B_coarse * controlPoints';
new_control_points(3:5,:) = new_control_points(1:3,:);
new_control_points(1,:) = controlPoints(:,1)';
new_control_points(2,:) = controlPoints(:,2)';

% make B-spline
pcrv = nrbmak((cat(2, new_control_points, complete_contolpoints))', target_knots);

% plot
nrbplot(pcrv, 200);

% refine back to the original mesh
% B_ref = basiskntins (p,target_knots,source_knots);
B_ref =  bzrproj_href( p, C, C_refined, 2);

% h-refinement
ref_control_points = B_ref * controlPoints';

% make B-spline
rcrv = nrbmak((cat(2, ref_control_points, [0 0 0 0 0 0 0 0 0 0]', [1 1 1 1 1 1 1 1 1 1]'))', refined_knots);

% plot
nrbplot(rcrv, 200);

%% Test Projector ---------------------------------------------------------

% Check if B_coares is the left inverse odf B_ref_test
B_ref_test =  bzrproj_href(p, C_coarse, C, 2);
Id = B_coarse*B_ref_test;

if (~isequal(normest(Id), 1))
    disp('is not a projection');
end

ref_test_control_points = B_ref_test * new_control_points;

% make B-spline
tcrv = nrbmak((cat(2, ref_test_control_points, [0 0 0 0 0 0]', [1 1 1 1 1 1]'))', source_knots);

% project control values of the source mesh onto the target mesh
B_coarse_test =  bzrproj_hcoarse( p, C_refined, C);

coarse_test_control_points = B_coarse_test * ref_control_points;

% make B-spline
ptcrv = nrbmak((cat(2, coarse_test_control_points, [0 0 0 0 0 0]', [1 1 1 1 1 1]'))', source_knots);

% plot
nrbplot(ptcrv, 200);
