function [epsilon_f, epsilon_r,capacity_f,capacity_r] = makeEpsilonSeq(model, reactions, epsilon0, k)
%default is take epsilon_f || epsilob_r = min(epsilon, K*Vmax_f || -K*Vmax_r)
% 06262019: the capacity output is added. this represent the flux carrying
%           capacity for each direction of a reaction
epsilon_f = epsilon0 * ones(length(reactions),1);
epsilon_r = epsilon0 * ones(length(reactions),1);
capacity_f = zeros(length(reactions),1);
capacity_r = zeros(length(reactions),1);

tol = 1e-7; % the default flux numerical tolerance for capacity output
for i = 1: length(reactions)
    testModel = model;
    testModel = changeObjective(testModel,reactions(i));
    Vm_f = optimizeCbModel(testModel, 'max');
    if Vm_f.obj*k < epsilon0 && Vm_f.obj*k > 1e-10 %avoid numeric error
        epsilon_f(i) = Vm_f.obj * k;
    end
    if Vm_f.obj > tol
        capacity_f(i) = true;
    else
        capacity_f(i) = false;
    end

    Vm_r = optimizeCbModel(testModel,'min');
    if -Vm_r.obj*k < epsilon0 && Vm_r.obj*k < -1e-10
        epsilon_r(i) = -Vm_r.obj * k;
    end
    if Vm_r.obj < -tol
        capacity_r(i) = true;
    else
        capacity_r(i) = false;
    end
end
