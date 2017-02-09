function [ w ] = wgheval_LLSQ( p, ndofs, index_dof )
%WGHEVAL_LLSQ: Bézier projection weights valuation
% Evaluate the weights to smooth the projected control values.
%
% Calling Sequence:
%
%   G_inv = grmninv( p, U, n )
%
%    INPUT:
%
%      p                    - polynomial degree
%      ndof                 - number of degrees of freedom
%      index_dof            - index of actual degree of freedom
%
%    OUTPUT:
%
%      w - weights of Bézier projection operator
%
%   Adapted from "Convergence of an
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

if (index_dof < p)
    w = 1/index_dof;
elseif (index_dof > ndofs - (p+1))
    w = 1/(ndofs-index_dof+1);
else
    w = 1/(p+1);
end

end

%! test: 
%! if (isequal(wgheval_LLSQ( 2, 5, 5 ), 1) && ...
%!      isequal(wgheval_LLSQ( 4, 10, 5 ), 1/5) && ...
%!       isequal(wgheval_LLSQ( 2, 4, 3 ), 1/2) && ...
%!        isequal(wgheval_LLSQ( 3, 6, 4 ),1/3))
%!     disp('Sucess !!!');
%! else
%!     disp('Error !!!');
%! end

