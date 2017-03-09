function f = moving_thermal_heat_source


f =@(x,y,xpath,ypath) 5830.0/(0.01*pi*0.015)...
    * exp(-3*(x-xpath).^2/0.015^2 - 3*(y-ypath).^2/0.01^2);


end

