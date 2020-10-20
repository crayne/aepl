%h = figure;
% hp = uipanel(...
%     'Title','Main Panel',...
%     'FontSize',12,...
%     'BackgroundColor','white',...
%     'Position',[.25 .1 .67 .67]...
% );
%          
% hsp = uipanel(...
%     'Parent',hp,...
%     'Title','Subpanel',...
%     'FontSize',12,...
%     'Position',[.4 .1 .5 .5]...
% );
% hbsp = uicontrol('Parent',hsp,'String','Push here',...
%               'Position',[18 18 72 36]);

%f = figure;
%f = gcf;
%clf

exper.path ='';
exper.plate_map = [];

close all
f = figure;
h = uitable();
h.Units = 'normalized';
h.Position = [.05,.05,.5,.9];

data = cell(20,2);
h.Data = data;
%h.RowName = 'butts';
h.RowName = fieldnames(exper);
h.ColumnEditable = true;


        %'CellEditCallback',@converttonum);

    h.CellEditCallback = @cellEdited;
    

h.Data(1,2) ={'<html><table border=0 width=400 bgcolor=#FF0000 color=#FF0000><TR><TD>Hello</TD></TR> </table></html>' }
h.Title('User provided')
function cellEdited(hObject,callbackdata)
     %numval = eval(callbackdata.EditData);
     r = callbackdata.Indices(1)
     c = callbackdata.Indices(2)
     newdata = callbackdata.EditData
     %hObject.Data{r,c} = numval; 
end







