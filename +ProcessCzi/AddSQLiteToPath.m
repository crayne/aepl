function AddSQLiteToPath()
    %CHECKJARPATH Summary of this function goes here
    %   Detailed explanation goes here

    %% ensure that the sqlite jar file is on the path
    dynamicPaths = javaclasspath('-dynamic');
    bfIsLoaded = false;
    if (~isempty(dynamicPaths))
        for i=1:length(dynamicPaths)
            [~,name,~] = fileparts(dynamicPaths{i});
            if (strcmpi('sqlite-jdbc-3.21.0',name))
                bfIsLoaded = true;
                break
            end
        end
    end

    if (~bfIsLoaded)
        curPath = mfilename('fullpath');
        [pathstr,~,~] = fileparts(curPath);
        javaaddpath(fullfile(pathstr,'sqlite-jdbc-3.21.0.jar'),'-end');
    end
end