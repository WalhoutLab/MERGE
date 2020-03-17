function cellWrite(filename,origCell)
% save a new version of the cell for reference
modCell = origCell;
% assume some cells are numeric, in which case set to char
iNum = cellfun(@isnumeric,origCell);
%% Method 1 using indexing and straight conversion = 7 sec
tic
modCell(iNum) = cellstr(num2str([origCell{iNum}]'));
% toc
%% Method 2 using arrayfun = 25 sec
% tic
% modCell(iNum) = arrayfun(@num2str,[origCell{iNum}],'unif',0);
% toc
% tic
[rNum,cNum] = size(origCell);
frmt = repmat([repmat('%s,',1,cNum-1),'%s\n'],1,rNum);
modCell = modCell';
fid = fopen(filename,'w');
fprintf(fid,frmt,modCell{:});
fclose(fid);
toc
end