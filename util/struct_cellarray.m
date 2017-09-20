function mycell = struct_cellarray(myds)

% handles conversion of a struct to a cell array when some of the struct's
% fields are numeric arrays

f = fieldnames(myds);
for k = 1:numel(f)
    if and(~iscell(myds.(f{k})), or(isnumeric(myds.(f{k})), islogical(myds.(f{k}))))
        myds.(f{k}) = num2cell(myds.(f{k}));
    end
end

mycell = struct2cell(myds);

if and(size(mycell,2) == 1, numel(myds.(f{1})) > 1)
    mycell = horzcat(mycell{:});
end
end