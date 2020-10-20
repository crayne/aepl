function [conditions, conditDict] = MakeConditDict(plateMapFile)

ttable = readtable(plateMapFile,'Delimiter',',','ReadVariableNames',false);


startRowIndices = find(strcmp(ttable{:,1}, 'B')); %T.Var1,'B'))

groupIndex = find(strcmp(ttable{:,1},'Groups')) + 2;
if ~isempty(groupIndex) 
%temp = grouping+2;
    %g = ;
    startRowIndices(find(startRowIndices(:) == groupIndex)) = [];
else
    groupIndex = false;
end

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


wellDict = containers.Map('KeyType','char','ValueType','char');
wellDict2 = containers.Map('KeyType','char','ValueType','char');


conditions = cell(0,2);
groups = cell(0,2);
%condits = Condit.empty();



tableCount = length(startRowIndices);
for r = 0:stopRowIndices(1)-startRowIndices(1)         
    for col = startColIndex:stopColIndex
        
        
        boop = cell(tableCount,1);
        for i = 1:tableCount
            boop{i} = ttable{startRowIndices(i)+r,col}{1};
        end
        if ~any(cellfun(@isempty,boop))     %if no elements of boop are empty
            boop2 = strjoin(boop,', ');
            wellDict(strcat(abc(r+1),oneTwoThree{col-1})) = boop2;
            

            if groupIndex
                wellDict2(strcat(abc(r+1),oneTwoThree{col-1})) = ttable{groupIndex+r,col}{1};
            end
        end
          
    end
end
%end
% 
% if length(startRowIndices) == 1 + length(grouping)
%     for row = startRowIndices(1):stopRowIndices(1)
%         for col = startColIndex:stopColIndex 
%             temp = ttable{row,col};
%             if ~(isempty(temp))
%                 wellDict(strcat(ttable{row,1}{1},oneTwoThree{col-1})) = temp;
%             end
%         end
%     end
% elseif length(startRowIndices) == 2 + length(grouping)
%     for r = 0:length(startRowIndices(1):stopRowIndices{1}) - 1
%         for col = startColIndex:stopColIndex 
%             
%             a = ttable{startRowIndices(1)+r,col}{1};
%             b = ttable{startRowIndices(2)+r,col}{1};
% %             disp(a)
% %             disp(isempty(a))
% %             disp(b)
% %             disp(isempty(b))
%             if ~(isempty(a) || isempty(b))
%                 temp = strcat(a,", ",b);
%                 
%             %disp(strcat('lsjfa;dlfj;ldas:',temp))
% 
%                 temp = char(temp(:));
%                 if ~(length(temp)==1)
%                     wellDict(strcat(abc(r+1),oneTwoThree{col-1})) = temp;
%                 end
%             end
%         end
%     end
% else 
%     disp("uh oh - too many tables?")
% end

    
    

conditions = unique(wellDict.values());

conditDict = containers.Map;%('KeyType','char','ValueType','char');

conditCount = size(conditions,2);
%condits = Condit.empty(conditCount,0);

for condit = conditions
    %disp(condit{1})
    %disp(class(condit{1}))
    conditDict(condit{1}) = {};
end
% 
% for i = 1:conditCount
%     condits(i).name = conditions(i);
% end

% for key = myMap2.keys()
%     %disp(key)
%     
% end
%wells = wellDict.keys();
%wells = mat2cell(wells);
for key = wellDict.keys() 
    %disp(class(key))
    %disp(myMap(key))
    well = key{1};
    condit = wellDict(well);
    %disp(condit)
    
    if conditDict.isKey(condit)
        %disp('good')
        len = length(conditDict(condit));
        %conditDict(condit) = {conditDict(condit),well};
        temp = conditDict(condit);
        %disp(temp)
        temp{len+1} = well;
        conditDict(condit) = temp;
        %disp(conditDict(condit))
        %disp(length(conditDict(wellDict(key{1}))))
        % = myMap(key{1})
    %if iskey(myMap2,myMap(key{1}))
%     else
%        disp(key) 
    end
end

% for key = conditDict.keys()
%     disp(key{1})
%     disp(conditDict(key{1}))
% end


end




