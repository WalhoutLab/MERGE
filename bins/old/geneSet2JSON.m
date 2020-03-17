function myjson = geneSet2JSON(data)
myjson = jsonencode(data);
for i = 1:length(myjson) 
    if myjson(i) == '"'
        myjson(i) = "'";
    end
end
end