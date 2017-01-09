function f = moving_heat_source_3D


f =@(x,y,z,t) 58.3e+05/(0.1*pi*0.15) * exp(-3*(x).^2/0.1^2 - 3*(y).^2/0.15^2);



end
