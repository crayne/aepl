%% ReadPlateMap  -  beta / almost stable
%   takes path of plate map
%   reads it in and
%   returns conditWellMap
%       map of 
%       conditions (str) -> wells (str[])
%

% need to prevent groups from being added to condition names



function [conditions,conditWellMap, groupConditMap] = ReadPlateMap(plateMapFile)

    disp(plateMapFile)
    ttable = readtable(plateMapFile,'Delimiter',',','ReadVariableNames',false);

    startRowIndices = find(strcmp(ttable{:,1}, 'B'));

    groupIndexTemp = find(strcmp(ttable{:,1},'Groups'));% + 2;
    if ~isempty(groupIndexTemp) 
        %disp(StartRowIndices
        %startRowIndices(find(groupIndex - 5 <= startRowIndices(:) <= groupIndex + 5)) = [];
        groupIndex = startRowIndices((startRowIndices(:) <= groupIndexTemp + 5) & (startRowIndices(:) >= groupIndexTemp));
        if length(groupIndex) == 1
        
        %temp = (startRowIndices <= groupIndexTemp + 5) && (startRowIndices >= groupIndexTemp);
        %group
            startRowIndices(startRowIndices(:) == groupIndex) = [];
        else
            groupIndex = false;
        end
    else
        groupIndex = false;
    end
    
    %if GroupIndex

    stopRowIndices = double.empty(length(startRowIndices),0);

    for z = 1:length(startRowIndices)
        abc = 'BCDEFGHIJKLMNOPQ';
        for i = 1:length(abc)   
            if i +startRowIndices(z) - 1 >= height(ttable)
                stopRowIndices(z) = i+startRowIndices(z)-1;
                break
            elseif ~strcmp(ttable{i+startRowIndices(z)-1,1}{1},abc(i))
                stopRowIndices(z) = i+startRowIndices(z)-2;
                break
            end

        end
    end

    
    oneTwoThree = {'2','3','4','5','6','7','8','9','10','11','12',...
        '13','14','15','16','17','18','19','20','21','22','23','24'};
    oneTwoThreeB = {'02','03','04','05','06','07','08','09','10','11','12',...
        '13','14','15','16','17','18','19','20','21','22','23','24'};

    startColIndex = 2;
    stopColIndex = size(ttable,2);
    for i = 2:size(ttable,2)    
        if ~strcmp(ttable{startRowIndices(z)-1,i},oneTwoThree{i-1})
            stopColIndex= i-1;
            break
        end
    end


    wellConditMap = containers.Map('KeyType','char','ValueType','char');
    wellGroupMap = containers.Map('KeyType','char','ValueType','char');

    tableCount = length(startRowIndices);
    for r = 0:stopRowIndices(1)-startRowIndices(1)         
        for col = startColIndex:stopColIndex

            temp = cell(tableCount,1);
            for i = 1:tableCount
                temp{i} = ttable{startRowIndices(i)+r,col}{1};
            end
            if ~any(cellfun(@isempty,temp))     %if no elements of temp are empty
                temp2 = strjoin(temp,', ');
                wellConditMap(strcat(abc(r+1),oneTwoThreeB{col-1})) = temp2;

                if groupIndex
                    wellGroupMap(strcat(abc(r+1),oneTwoThreeB{col-1})) = ttable{groupIndex+r,col}{1};
                end
            end
        end
    end


    conditions = unique(wellConditMap.values());
    groups = unique(wellGroupMap.values());

    conditWellMap = containers.Map;
    for condit = conditions
        conditWellMap(condit{1}) = {};
    end

    groupConditMap = containers.Map;
    for group = groups
        groupConditMap(group{1}) = {};
    end


    for key = wellConditMap.keys()
        well = key{1};
        condit = wellConditMap(well);
        conditWellMap(condit) = [conditWellMap(condit), well];

        try 
            group = wellGroupMap(well);
            if ~any(strcmp(groupConditMap(group),condit))
                groupConditMap(group) = [groupConditMap(group),condit];
            end
        catch             
        end
    end

end




