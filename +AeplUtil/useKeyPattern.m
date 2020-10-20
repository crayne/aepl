


function value = useKeyPattern(map,pattern) 
    key = AeplUtil.FindKey(map,pattern);
    if iscell(key) && length(key) == 1 
        key = key{1};
    end
    if ~key
        error('The specified key is not present in this container.')
    end
    value = map(key);
    if iscell(value) && length(value) == 1 
        value = value{1};
    end
       

end