function expression_logic = freezeANDlogic(expression)
NONETYPE = '0';% by default, the undetected genes are NaN (will be ignored)
NUMBER = '[0-9\.\-e]+';
MAYBE_NUMBER = [NUMBER '|' NONETYPE];
%freeze the AND logic genes
expression = regexprep(expression,'and','&');
expression = regexprep(expression,'or','|');
numStarts = regexp(expression,MAYBE_NUMBER);
numEnds = regexp(expression,MAYBE_NUMBER,'end');
if numStarts(1) > 1
    expression_logic = expression(1:numStarts(1)-1);
else
    expression_logic = '';
end
for i = 1:length(numStarts)-1
    %replace and judge 
    tmp = [expression(1:numStarts(i)-1),'0',expression(numEnds(i)+1:end)];
    if ~eval(tmp) %if this is "AND" logic gene
        expression_logic = [expression_logic, 'fz_',expression(numStarts(i):numEnds(i)),'_fz',expression(numEnds(i)+1:numStarts(i+1)-1)];
    else
        expression_logic = [expression_logic, expression(numStarts(i):numStarts(i+1)-1)];
    end
end
%the last number
if numEnds(end) < length(expression)
    tmp = [expression(1:numStarts(end)-1),'0',expression(numEnds(end)+1:end)];
    if ~eval(tmp) %if this is "AND" logic gene
        expression_logic = [expression_logic, 'fz_', expression(numStarts(end):numEnds(end)),'_fz',expression(numEnds(end)+1:end)];
    else
        expression_logic = [expression_logic, expression(numStarts(end):end)];
    end
else
    tmp = [expression(1:numStarts(end)-1),'0'];
    if ~eval(tmp) %if this is "AND" logic gene
        expression_logic = [expression_logic, 'fz_', expression(numStarts(end):numEnds(end)),'_fz'];
    else
        expression_logic = [expression_logic, expression(numStarts(end):end)];
    end
end
expression_logic = regexprep(expression_logic,'&','and');
expression_logic = regexprep(expression_logic,'\|','or');