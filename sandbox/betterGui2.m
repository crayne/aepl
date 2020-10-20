


close all
gui.f = figure;

% <sizeSettings>

    % si = sizes
    si.s = .03;      % s = spaces (padding)

    si.p1h = 3/5;    % p1h = panel1 height
    si.p1w = 1/3;    % p1w = panel1 width

    si.pch = 4/5;    % pch = processCziPanel height

    si.sdh = .6;     % sdh = summarizeDataPanel height
    si.sdw = 2/5;    % sdw = summarizeDataPanel width



    si.gpph = .06;      % gpph = getPathPanel height
    si.gptw = 1/8;      % gptw = .getPathText width

    si.fs = 12;         % fs = fontSize
    
% </sizeSettings>
% <getPath>

    % gpp = getPathPanel
    gui.gpp = uipanel(gui.f, 'Position',[si.s, 1-si.s-si.gpph, 1-2*si.s, si.gpph]);

        % gpt = getPathText
        gui.gpt = uicontrol(gui.gpp, 'Style','pushbutton', 'Units','normalized',...
            'FontSize',si.fs, 'HorizontalAlignment','left', 'String','Path:',...
            'Callback',@gptCallback);
            gui.gpt.Position = [0, si.s, si.gptw, gui.gpt.Position(4)];

        % gpe = getPathEdit
        gui.gpe = uicontrol(gui.gpp, 'Style','edit', 'Units','normalized',...
            'FontSize',si.fs, 'HorizontalAlignment','left',...
            'Position',[gui.gpt.Position(3), si.s, 1-gui.gpt.Position(3), gui.gpt.Position(4)],...
            'Callback',@gpeCallback);

% </getPath>
% <table1>

    % t1 = table1
    gui.t1 = uitable(gui.f, 'Units','normalized', 'FontSize',si.fs,...
        'Position',[si.s, si.p1h+si.s, 1-2*si.s, 1-si.p1h-3*si.s-si.gpph]);

% </table1>
% <panel1>

    % p1 = panel1
    gui.p1 = uipanel(gui.f, 'Position',[si.s, si.s, si.p1w-si.s, si.p1h-si.s]); 
        
        % czir = czi required
        gui.czir = uicontrol(gui.p1, 'Style','text', 'Units','normalized',...
            'FontSize',si.fs, 'HorizontalAlignment','left', 'String','Requires czi files',...
            'Position',[2*si.s, si.pch+si.s, 1-3*si.s, 1-2*si.s-si.pch]);
         
        % pc = ProcessCzi (Panel)
        gui.pc = uipanel(gui.p1, 'FontSize',si.fs,...
            'Position',[si.s, si.s, 1-2*si.s, si.pch-si.s], 'Title','ProcessCzi');
            % pcRun = ProcessCziRun (Button)
            gui.pcRun = uicontrol(gui.pc, 'Units','normalized', 'Style','pushbutton',...
                'FontSize',si.fs, 'String','Run');
                gui.pcRun.Position = [si.s,si.s,gui.pcRun.Position(3),gui.pcRun.Position(4)];
            
            % bg = buttongroup
            gui.bg = uibuttongroup(gui.pc, 'Units','normalized',...
                'Position',[si.s, gui.pcRun.Position(4)+2*si.s, 1-2*si.s, 1-gui.pcRun.Position(4)-3*si.s]);

                gui.notif = uicontrol(gui.bg, 'Style','radiobutton', 'Units','normalized',...
                    'FontSize',si.fs, 'String', '<html>Create tracking<br>tifs</html>');
                si.rbh = gui.notif.Position(4)*2;
                gui.notif.Position = [si.s, 1-si.s-si.rbh ,1-2*si.s,si.rbh];
                gui.tif = uicontrol(gui.bg, 'Style','radiobutton', 'Units','normalized',...
                    'FontSize',si.fs,'Position',[si.s, 1-2*si.s-2*si.rbh,1-2*si.s,si.rbh],...
                    'String','<html>Do not create<br>tracking tifs</html>');
% </panel1>
% <panel2>

    % p2 = panel2
    gui.p2 = uipanel(gui.f, 'Position',[si.p1w+si.s, si.s, 1-2*si.s-si.p1w, si.p1h-si.s]);

        % csvr = csv required
        gui.csvr = uicontrol(gui.p2, 'Style','text', 'Units','normalized',...
            'FontSize',si.fs, 'HorizontalAlignment','left','String','Requires csv files');
            si.r2h = gui.csvr.Position(4);
            gui.csvr.Position = [si.s, 1-2*si.s-si.r2h, si.sdw-2*si.s, si.r2h]; 

        % pmr = platemap required
        gui.pmr = uicontrol(gui.p2, 'Style','text', 'Units','normalized', 'FontSize',si.fs,...
            'HorizontalAlignment','left', 'String','Requires platemap file',...
            'Position',[si.s, 1-3*si.s-2*si.r2h, si.sdw-2*si.s, si.r2h]);



        % sd = summarizeData (Panel)
        gui.sd = uipanel(gui.p2, 'FontSize',si.fs, 'Title','SummarizeData',...
            'Position',[si.s, si.s, si.sdw-2*si.s, si.sdh-si.s]);
            
            % sdRun = summarizeDataRun (Button)
            gui.sdRun = uicontrol(gui.sd, 'Style','pushbutton', 'Units','normalized',...
                'FontSize',si.fs, 'String','Run');
            gui.sdRun.Position = [si.s,si.s,gui.sdRun.Position(3),gui.sdRun.Position(4)];

        % pd = plotData (Panel)
        gui.pd = uipanel(gui.p2, 'FontSize',si.fs, 'Title','PlotData',...
            'Position',[si.sdw, si.s, 1-si.sdw-si.s, 1-2*si.s]);
        
            % cSet = byConditionSetLayout
            gui.cSet = uicontrol(gui.pd, 'Style','checkbox', 'Units','normalized',...
                'FontSize',si.fs, 'String','byConditionSetLayout');
            si.cbh = gui.cSet.Position(4)*2;
            gui.cSet.Position = [si.s, 1-si.s-si.cbh, 1-2*si.s, si.cbh];
            
            % pdRun = plotDataRun (Button)
            gui.pdRun = uicontrol(gui.pd, 'Style','pushbutton', 'Units','normalized',...
                'FontSize',si.fs, 'String','Run');
                gui.pdRun.Position = [si.s,si.s,gui.pdRun.Position(3),gui.pdRun.Position(4)];

            % grp = groupsRequiredPanel
            gui.grp = uipanel(gui.pd);
            gui.grp.Position = [si.s, gui.pdRun.Position(4)+2*si.s , 1-2*si.s, 1-3*si.s-si.cbh- gui.pdRun.Position(4)];

                % gr = groupsRequired (text)
                gui.gr = uicontrol(gui.grp, 'Style','text', 'Units','normalized', 'FontSize',si.fs,...
                    'HorizontalAlignment','left',  'String','Requires groups in platemap file');


                gui.gr.Position = [si.s, 1-si.s-2*gui.gr.Position(4), 1-2*si.s, 2*gui.gr.Position(4)];

                % cAuto = byConditionAutoLayout
                gui.cAuto = uicontrol(gui.grp, 'Style','checkbox', 'Units','normalized',...
                    'FontSize',si.fs, 'String','byConditionAutoLayout',...
                    'Position',[2*si.s, 1-gui.gr.Position(4)-2*si.s-si.cbh, 1-3*si.s, si.cbh]);
                % g = byGroup
                gui.g = uicontrol(gui.grp, 'Style','checkbox', 'Units','normalized',...
                    'FontSize',si.fs, 'String','byGroup',...
                    'Position',[2*si.s, 1-gui.gr.Position(4)-3*si.s-2*si.cbh, 1-3*si.s, si.cbh]);

           

% </panel2>
% <functions>



    function gptCallback(hObject,callbackdata)
        
        startin = '/Volumes/baylieslab/Current Lab Members/Whitney/Rhabdomyosarcoma plate movies/';
        temp = uigetdir(startin);
        if ~strcmp(temp,0)
            hObject.Parent.Children(1).String = temp;
            gpeCallback(hObject.Parent.Children(1),[])
        end
    end



    function gpeCallback(hObject,callbackdata)
         
        hasCzi = checkForCzi(hObject.String);
       
        if startsWith(hasCzi, 'czi found:')
            t = split(hasCzi,':');
            f = hObject.Parent.Parent;
           
            
            pcp = f.Children(2).Children(1);
            
            
            gree = [0.2,0.55,0.1];
            
            czir =  f.Children(2).Children(2);
            czir.ForegroundColor = gree;
           
            czir.String = [czir.String;strcat("  found ",t(2), " czi files")];
            
            pcp.Children(2).ForegroundColor = gree;

        end

    end


    function hasCzi = checkForCzi(experPath)
        if ~isfolder(experPath)
            hasCzi = 'cannot find folder';
        
        else
        
            contents = dir(experPath);
            count = 0;
            for i = 1:length(contents)
                if endsWith(contents(i).name,'.czi') 
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
    
    
% </functions>
    
    
    

