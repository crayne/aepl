%% FindFile  -  dev
%   helper function
%       takes folder and string pattern
%       returns first file in folder containing pattern



function fullpath = FindFile(folder, pattern)
   

    contents = dir(folder);
    for i = 1:length(contents)
        if contains(contents(i).name,pattern)
           fullpath = fullfile(folder,contents(i).name);
           return
        end
    end
    try 
        error('Could not find file:\n\tfolder = %s\n\tpattern = %s',folder,pattern)

    catch e
        disp(folder)
        disp(pattern)
        error(e.getReport())

    end
    
    
end