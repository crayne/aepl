


function key = FindKey(map, pattern) 
    for k = map.keys() 
        if contains(k,pattern)
            key = k;
            return
        end
        
    end

    key = false;

end