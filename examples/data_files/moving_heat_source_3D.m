function f = moving_heat_source_3D

f = @(x,y,z,path_x,path_y,path_z) 6.0*sqrt(3.0)*58.3e+02/(pi * sqrt(pi) * 0.015 * 0.010 * 0.002) * ...
    exp(- 3*(x-path_x).^2/0.015^2 - 3*(y-path_y).^2/0.010^2 - 3*(z-path_z).^2/0.002^2);

end


