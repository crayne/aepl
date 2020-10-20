PLATE_MAP_SUF = '_plate-map.csv';
CZI_SUF = '.czi';


%% Add MATLAB Utilities for BioFormats 
addpath('../Utilities/src/MATLAB')

%addpath('../aeplClasses')


%% Directory 
curDir = pwd;
guiGetDir = @uigetdir;

f18_03_18 = '/Users/baylieslab/Documents/Amelia/rmsMim/18-03-18/';
f18_06_20 = '/Volumes/baylieslab/Current Lab Members/Whitney/Rhabdomyosarcoma plate movies/18-06-20/';


f16_06_23 = '/Volumes/baylieslab/Current Lab Members/Whitney/Rhabdomyosarcoma plate movies/16-06-23 (1st plate)/16-06-23/';

exper = struct();
%exper.folder = f16_06_23;

exper.folder = guiGetDir('/Volumes/baylieslab/Current Lab Members/Whitney/Rhabdomyosarcoma plate movies/');

%%
tic

[~,exper.name] = fileparts(exper.folder(1:end-1));


exper.plateMapFile = fullfile(exper.folder,[exper.name,PLATE_MAP_SUF]);

[exper.conditions,exper.conditDict] = MakeConditDict(exper.plateMapFile);


exper.czi = dir(fullfile(exper.folder,'*'+CZI_SUF));
exper.figureHandles = []; % Necessary?
exper.frameCount = 0;
% exper.wells = containers.Map('KeyType','char','ValueType','struct');



for condition = exper.conditions
    disp(condition{1})
    
    wells = exper.conditDict(condition{1});
    conditVeloc = cell(size(wells,2));
    buttts = true;
    for w = wells 
        
        well = struct();
        well.condition = condition{1};
        well.name = w{1};
        well.filename = strcat(exper.name,'_',well.name,CZI_SUF);
        well.path = strcat(exper.folder,well.filename);
            
        
            [P,F,~] = fileparts(well.path);
            P2 = fullfile(P,'Xls');
            if ~exist(P2,'dir')
                mkdir(P2)
            end
            XlsName = fullfile(P2,[F,'.csv']);
            
            
        if exist(well.path, 'file') %&& ~exist(XlsName,'file')
            try
                well = RunOneWell(exper,well);

            catch e
                fprintf(2,"condition: " + well.condition + ", well:" + well.name)
                fprintf(2,"exception: " + getReport(e)+"\n")
                
%                 disp("except.message: " + e.message)
%                 disp(strcat('exception.stack.line: ',e.stack.file))
                continue
                
                
    
            end
            

        end
    
    end

toc
end 


%% well = struct with fields
%   condition
%   filename
function well = RunOneWell(exper,well)
                disp(well.name)

    [well.im,well.imd] =  MicroscopeData.Original.ReadData(exper.folder,well.filename);

    well.imdim = size(well.im);
    if exper.frameCount == 0 
        exper.frameCount = well.imdim(5);
    end

    [well.cells] = segIms(well.im);


    [well.cells2,well.edges] = getTracks(well.cells,well.imdim);

    
    %[P,F,~] = fileparts(Name);
    P2 = fullfile(exper.folder,'outTifs');
    if ~exist(P2,'dir')
        mkdir(P2)
    end
    %XlsName = fullfile(P2,[F,'.csv']);
    well.outTifPath = fullfile(P2,strcat(exper.name,'_',well.name,'.tif'));
    DrawTracks(squeeze(well.im),well.cells2,well.outTifPath);
    
    ExportTrackStats(well.cells2,well.imdim,well.path)
            
    well.allSegs = vertcat(well.cells2{:});
    well.allTracks = [well.allSegs.Tid];
    well.tracks = unique(well.allTracks);
    well.cellCount = size(well.tracks,2);
    disp("well.cellCount = " + well.cellCount)
    well.velocs = cell(exper.frameCount,size(well.tracks,2));

    %ShowPlot(well.cells2,well.name,well.condition,exper.figureHandles)
    
end


