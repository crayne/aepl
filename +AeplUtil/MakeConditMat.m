%% MakeConditMat  -  alpha
%   takes a condit and gets all displacement columns
%   returns condit with new field condit.mat - matrix of displacement cols
%
%
%



function condit = MakeConditMat(condit)
    condit.mat = [];
    col = 1;
    for w = 1:length(condit.wells)
        try 
            %disp(w)
            if condit.wells(w).path == 0
                continue
            end
                
            
            frames = length(condit.wells(w).cells(1).distance(2:end));

            condit.wells(w).mat = [];
            wCol = 1;
            
            for c = 1:length(condit.wells(w).cells)
                try 
                    condit.mat(1:frames,col) = condit.wells(w).cells(c).distance(2:end);
                    col = col+1;


                    condit.wells(w).mat(1:frames,wCol) = condit.wells(w).cells(c).distance(2:end);
                    wCol = wCol + 1;

                    condit.wells(w).cells(c).coords = horzcat(condit.wells(w).cells(c).xcoords, condit.wells(w).cells(c).ycoords);


                catch e
                    fprintf(2,['condition: ', condit.name,'\n'])
                    fprintf(2,[e.getReport(),'\n'])
                end
            end
        catch e
            fprintf(2,['condition: ', condit.name,'\n'])
            fprintf(2,[e.getReport(),'\n'])
        end
        
        
        
        
    end
    condit.mat(condit.mat>50)=nan;
end
