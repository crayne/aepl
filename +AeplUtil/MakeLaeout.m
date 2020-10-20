%% MakeLaeout  -  dev
%   given number of plots returns a layout 
%       (laeout for indentifier uniqueness)
%   from montage maker code in fiji
%   gives row and column counts close-ish to a square



function laeout = MakeLaeout(plotCount)
    if plotCount == 1
       laeout = [3,4]; 
   
    else
        colCount = floor(sqrt(plotCount));
        rowCount = colCount;
        n = plotCount - colCount*rowCount;
        if n > 0 
            colCount = colCount + ceil(n/rowCount);
        end
        %laeout = [rowCount colCount];
        laeout = [colCount,rowCount];
    end
end
