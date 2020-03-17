function myjson = epsilon2JSON(model,epsilon_f,epsilon_r)
data = struct();
for i = 1:length(model.rxns)
    data.(model.rxns{i}) = [epsilon_f(i),epsilon_r(i)];
end
myjson = jsonencode(data);
for i = 1:length(myjson) 
    if myjson(i) == '"'
        myjson(i) = "'";
    end
end
end