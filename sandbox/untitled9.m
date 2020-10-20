


close all
f = figure;

s = .03;

p1h = 3/5;

p1w = 1/3;

pch = 4/5;

sdh=.6;
sdw = 1/2;



lh = .05;
fs = 12;

p1 = uipanel(f, 'Position',[s, s, p1w-s, p1h-s]);
p2 = uipanel(f, 'Position',[p1w+s, s, 1-2*s-p1w, p1h-s]);

t1 = uitable(f, 'Units','normalized', 'FontSize',fs, 'Position',[s, p1h+s, 1-2*s, 1-p1h-3*s-lh]);

gpw = 1/7;
gpp = uipanel(f, 'Position',[s, 1-s-lh, 1-2*s, lh]);
gpt = uicontrol(gpp, 'Style','text', 'Units','normalized', 'FontSize',fs, 'Position',[0, 0, gpw, 1], 'String','Path:');
getPath = uicontrol(gpp, 'Style','edit', 'Units','normalized', 'FontSize',fs, 'Position',[gpw, 0, 1-gpw, 1]);



pc = uipanel(p1, 'FontSize',fs, 'Position',[s, s, 1-2*s, pch-s], 'Title','ProcessCzi');
%cw = uicontrol(p1,'Style','text', 'Position',[s, pch+s, 1-s, 1-3*s-pch], 'String','Requires czi files');
%cwp = uipanel(p1,'Position',[s, pch+s, 1-2*s, 1-2*s-pch]);
czir = uicontrol(p1, 'Style','text', 'Units','normalized', 'FontSize',fs, 'HorizontalAlignment','left', 'Position',[2*s, pch+s, 1-3*s, 1-2*s-pch] ,'String','Requires czi files');



sd = uipanel(p2, 'FontSize',fs, 'Position',[s, s, sdw-2*s, sdh-s], 'Title','SummarizeData');

pd = uipanel(p2, 'FontSize',fs, 'Position',[sdw, s, sdw-s, 1-2*s], 'Title','PlotData');


%r2h = (1-sdh-3*s)/2;
csvr = uicontrol(p2, 'Style','text', 'Units','normalized', 'FontSize',fs,...
    'HorizontalAlignment','left','String','Requires csv files');
r2h = csvr.Position(4);
csvr.Position = [s, 1-2*s-r2h, sdw-2*s, r2h]; 

pmr = uicontrol(p2, 'Style','text', 'Units','normalized', 'FontSize',fs,...
    'HorizontalAlignment','left', 'Position',[s, 1-3*s-2*r2h, sdw-2*s, r2h] ,'String','Requires platemap file');


sdRun = uicontrol(sd, 'Style','pushbutton', 'Units','normalized', 'FontSize',fs, 'String','Run');
sdRun.Position = [s,s,sdRun.Position(3),sdRun.Position(4)];


pdRun = uicontrol(pd, 'Style','pushbutton', 'Units','normalized', 'FontSize',fs, 'String','Run');
pdRun.Position = [s,s,pdRun.Position(3),pdRun.Position(4)];



%pcRun = uicontrol(pc, 'Units','normalized', 'Style','pushbutton', 'FontSize',fs, 'Position',[s,s,1-2*s,.3], 'String','Run');
pcRun = uicontrol(pc, 'Units','normalized', 'Style','pushbutton', 'FontSize',fs, 'String','Run');
%p = pcRun.Position;
pcRun.Position = [s,s,pcRun.Position(3),pcRun.Position(4)];

bg = uibuttongroup(pc, 'Units','normalized', 'Position',[s,pcRun.Position(4)+2*s,1-2*s,1-pcRun.Position(4)-3*s]);

notif = uicontrol(bg, 'Style','radiobutton', 'Units','normalized', 'FontSize',fs, 'String','Create tracking tifs');
rbh = notif.Position(4)*2;
notif.Position = [s, 1-s-rbh ,1-2*s,rbh];
tif = uicontrol(bg, 'Style','radiobutton', 'Units','normalized', 'FontSize',fs,'Position',[s, 1-2*s-2*rbh,1-2*s,rbh], 'String','<html>Do not create<br>tracking tifs</html>');


cSet = uicontrol(pd, 'Style','checkbox', 'Units','normalized', 'FontSize',fs, 'String','cSet');
cbh = cSet.Position(4)*2;
cSet.Position = [s, 1-s-cbh, 1-2*s, cbh];
% cSet.Position(2) = 1-s-cbh;
% cSet.Position(1) = s;


grp = uipanel(pd);
%grp.Position(4) = 1-2*s-cSet.Position(4);
grp.Position = [s, pdRun.Position(4)+2*s , 1-2*s, 1-3*s-cbh- pdRun.Position(4)];

gr = uicontrol(grp, 'Style','text', 'Units','normalized', 'FontSize',fs, 'HorizontalAlignment','left',  'String','Requires groups in platemap file');
% gr.Position(1) = s;
% gr.Position(2) = 1-s-gr.Position(4);

gr.Position = [s, 1-s-2*gr.Position(4), 1-2*s, 2*gr.Position(4)];

cAuto = uicontrol(grp, 'Style','checkbox', 'Units','normalized', 'FontSize',fs,...
    'Position',[2*s, 1-gr.Position(4)-2*s-cbh, 1-3*s, cbh], 'String','cAuto');
g = uicontrol(grp, 'Style','checkbox', 'Units','normalized', 'FontSize',fs,...
    'Position',[2*s, 1-gr.Position(4)-3*s-2*cbh, 1-3*s, cbh], 'String','g');




% [s, 1-s-h, , ]



