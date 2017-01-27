function C = capacity(eu, eu_old )
%CAPACITY for Ti-6Al-4V from K. C. Mills,"Recommended Values of
%Thermophysic Properties for selected commercial alloys", 2002.

C = zeros(size(eu));

for i=1:size(eu,1)
    for j=1:size(eu,2)
        capacity_pc = 7500 * 600 + ...
            7500 * 261e+03 * func_pc_der(eu(i,j), eu_old(i,j)) * (tanh( (eu(i,j) - 800) / 800) + 1);
        
        C(i,j) =  capacity_pc;
    end
end



end

function f_pcDer = func_pc_der(eu, eu_old)

if eu == eu_old
    f_pcDer = 0.0;
else
    f_pcDer = (func_pc(eu) - func_pc(eu_old)) / (eu - eu_old);
end

end

function f_pc = func_pc(eu)

f_pc = zeros(size(eu));

f_pc(eu <= 800) = 0.0;
f_pc(eu > 800) = 1.0;

end

