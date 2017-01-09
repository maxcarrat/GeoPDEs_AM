%ADAPTIVITY_ESTIMATE_GRADIENT: Gradient base error indicators for non-linear Poisson transient problem, using globally smooth (C^1) hierarchical spaces.
%
% We consider the diffusion problem
%
%   c(x,u) du/dt - div ( epsilon(x,u) grad (u)) = f    in Omega = F((0,1)^n)
%                            epsilon(x,u) du/dn = g    on Gamma_N
%                                             u = h    on Gamma_D
%
% USAGE:
%
%   est = adaptivity_estimate_gradient (u, hmsh, hspace, problem_data, adaptivity_data)
%
% INPUT:
%
%   u:            degrees of freedom
%   hmsh:         object representing the hierarchical mesh (see hierarchical_mesh)
%   hspace:       object representing the space of hierarchical splines (see hierarchical_space)
%   adaptivity_data: a structure with the data for the adaptivity method. In particular, it contains the fields:
%    - flag:          'elements' or 'functions', depending on the refinement strategy.
%    - C0_est:        multiplicative constant for the error indicators
%
%
% OUTPUT:
%
%   est: computed a gradient based indicators

function est = adaptivity_estimate_gradient (u, time_step, hmsh, hspace, problem_data, adaptivity_data)

if (isfield(adaptivity_data, 'C0_est'))
    C0_est = adaptivity_data.C0_est;
else
    C0_est = 1;
end

[ders, F] = hspace_eval_hmsh (u, hspace, hmsh, {'gradient'});
dernum = ders{1};

x = cell (hmsh.rdim, 1);
for idim = 1:hmsh.rdim;
    x{idim} = reshape (F(idim,:), [], hmsh.nel);
    x{idim} =  x{idim} - ones(size(x{idim})).*problem_data.path(time_step, idim);
end

valf = problem_data.f(x{:}, problem_data.time_discretization(time_step+1));
aux = reshape (sum (C0_est .* dernum, 1), size(valf)).^2;

w = [];
h = [];
for ilev = 1:hmsh.nlevels
    if (hmsh.msh_lev{ilev}.nel ~= 0)
        w = cat (2, w, hmsh.msh_lev{ilev}.quad_weights .* hmsh.msh_lev{ilev}.jacdet);
        h = cat (1, h, hmsh.msh_lev{ilev}.element_size(:));
    end
end
h = h * sqrt (hmsh.ndim);

est = sqrt (sum (aux.*w));
est = C0_est*h.*est(:);


end