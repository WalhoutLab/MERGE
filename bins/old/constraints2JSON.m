function myjson = constraints2JSON(model)
data = struct();
for i = 1:length(model.rxns)
    data.(model.rxns{i}) = [model.lb(i),model.ub(i)];
end
myjson = jsonencode(data);
for i = 1:length(myjson) 
    if myjson(i) == '"'
        myjson(i) = "'";
    end
end
end