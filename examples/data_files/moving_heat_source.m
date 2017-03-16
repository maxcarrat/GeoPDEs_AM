function f = moving_heat_source

f =@(x,y,path_x,path_y)  3.0 * 5830.0/(pi * 0.0015 * 0.0010 ) * ...
    exp(- 3*(x-path_x).^2/0.0015.^2 - 3*(y-path_y).^2/0.0010.^2 );

end
