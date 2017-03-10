function f = moving_thermal_heat_source


f =@(x,y,xpath,ypath) 3.0 * 5830.0/(0.015*pi*0.010)...
    * exp(-3.0*(x-xpath).^2/0.015.^2 - 3.0*(y-ypath).^2/0.010.^2);


end

