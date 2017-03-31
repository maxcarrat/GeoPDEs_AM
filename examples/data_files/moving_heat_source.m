function f = moving_heat_source

f =@(x,y,path_x,path_y)  3.0 * 3.0 * 58300.0 * 5/(pi * 0.00015 * 0.00010 ) * ...
    exp(- 3*(x-path_x).^2/0.00015.^2 - 3*(y-path_y).^2/0.00010.^2 );

end
