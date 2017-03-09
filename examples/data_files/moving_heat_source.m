function f = moving_heat_source

f =@(x,y,path_x,path_y)  3.0 * 583000.0/(pi * 0.015 * 0.010 ) * ...
    exp(- 3*(x-path_x).^2/0.015.^2 - 3*(y-path_y).^2/0.010.^2 );

end
