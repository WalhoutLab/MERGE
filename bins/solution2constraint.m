function [MILProblem] = solution2constraint(MILProblem,solution)
obj = solution.obj;
A = MILProblem.A;
%add constraint
A(end+1,:) = MILProblem.c;
MILProblem.A = A;
%update others
MILProblem.b = [MILProblem.b;obj];
if MILProblem.osense == -1
    mycsense = 'G';
else
    mycsense = 'L';
end
if size(MILProblem.csense,1) == 1 %for some model, the csense is a string instead of vector
    csense = [MILProblem.csense,mycsense];
else %is a vector
    csense = [MILProblem.csense; mycsense];
end
MILProblem.csense = csense;
end