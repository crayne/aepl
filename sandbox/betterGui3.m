betterGui2

function gptCallback(hObject,callbackdata)
    hObject.Parent.Children(1).String = uigetdir();
    gpeCallback(hObject.Parent.Children(1),[])
end

% function getPathCallback(hObject,callbackdata)
% 
% end

function gpeCallback(hObject,callbackdata)
    disp('butts')
    hasCzi = checkForCzi(hObject.String);
    h = hObject.Parent;
    h.Children(2).String = 'butts';

    if startsWith(hasCzi, 'czi found:')
        %h = hObject.Parent.Parent;
        
        
    end

end


function hasCzi = checkForCzi(experPath)
    %bool = true;
    if ~isfolder(experPath)
        hasCzi = 'cannot find folder';
        
    else
        contents = dir(experPath);
        count = 0;
        for c = contents
            if endsWith(c.name,'.czi') 
                count = count + 1;
            end
                
        end
        if count > 0 
            hasCzi = ['czi found:',num2str(count)];
        else 
            hasCzi = ['no czi found'];
        end
        
    end

end