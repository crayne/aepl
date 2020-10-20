%% RunProcessCzi  -  beta / almost stable
%     wrapper function that calls the heavy lifters for processing czi
%
%   Does not require plate map
%

function RunProcessCzi(experPath)
disp('Starting RunProcessCzi')

%% Add MATLAB Utilities for BioFormats SC Maybe this is not necessary
addpath('./Utilities/src/MATLAB')

if ~exist('experPath','var')
    experPath = uigetdir();
end

global CONST
if isempty(CONST)
    ProcessCzi.initConstants
end

%% intialize exper
exper = struct();


dlist = dir(fullfile(experPath,['*',CONST.CZI_SUF]));
dlist = [ dlist;dir(fullfile(experPath,['*.tif']))];

exper.czi = dlist;
%exper.frameCount = 0;

S = split(experPath,'*');
experPath = S{1};
%% Drawing And Visualization

DrawPlot = CONST.DRAW_MODE;

Draw6Panel = CONST.DRAW_MODE;

if ~(CONST.DRAW_MODE)
    f = waitbar(0,'Please wait...');
    v = [];
    v2 = [];
else

    v = VideoWriter(fullfile(experPath,'OutputVideo.mp4'),'MPEG-4');
    v.Quality = 100;
    v.FrameRate = 15;
    open(v)


    v2 = VideoWriter(fullfile(experPath,'OutputVideo6Panel.mp4'),'MPEG-4');
    v2.Quality = 100;
    v2.FrameRate = 15;
    open(v2)

end
%%


for w = 1:length(dlist)

    Path1 = dlist(w).folder;

    %% Make tif dir and csv dir if they don't exist
    tifSaveDir = fullfile(Path1,CONST.ANNOTATED_TIF_DIR);
    if ~exist(tifSaveDir,'dir')
        mkdir(tifSaveDir)
    end
    %AeplUtil.labelVersion(tifSaveDir)

    csvSaveDir = fullfile(Path1,CONST.CSV_DIR);
    if ~exist(csvSaveDir,'dir')
        mkdir(csvSaveDir)
    end
    %         try
    %% parse name
    wellFile = exper.czi(w).name;
    temp = strsplit(wellFile,'.');
    temp = strsplit(temp{1},'-');

    wellName = temp{end};

    %% initialize paths for saving well data
    wellTifSavePath = fullfile(tifSaveDir,strcat(wellName,'.tif'));
    wellCsvSavePath = fullfile(csvSaveDir,strcat(wellName,'.csv'));
    wellBW_SavePath = fullfile(tifSaveDir,strcat(wellName,'BW.tif'));
    % used when doing one frame and DrawTracks2
    % wellTifSavePath2 = fullfile(tifSaveDir,strcat(wellName,'_.tif'));
    if ~(CONST.DRAW_MODE)
        waitbar(w/length(dlist),f,['Processing ',wellName]);
    end

    %% Process well czi
    %   only if a csv for that well doesn't already exist
    if exist(wellCsvSavePath,'file')
        disp(['Skipping well ', wellName, ', csv file already exists']);
        %continue;
    end

    disp(['Starting well ',wellName]);
    %% Open and load Image
    tic
    [im,imd] =  MicroscopeData.Original.ReadData(Path1,wellFile);
    imdim = size(im);
    fprintf(1,'\t\t')
    toc

    %% The heavy lifting - segment (or process) the image (series)
    tic
    [cells,BW] = ProcessCzi.SegIms(im,Draw6Panel,v2);

    %% Write Binaries to a file
    if 1 && ~isempty(BW) && DrawPlot
        delete(wellBW_SavePath)
    for i = 1:size(BW,3)

        imwrite(BW(:,:,i),wellBW_SavePath,'writemode','append')
    end
    end



    fprintf(1,'\t\t')
    toc
    %% from the raw data about every segmented cell from every frame,
    %   get tracks for each cell over all frames
    if isempty(cells)
        continue
    end

    tic
    [cells2,edges] = ProcessCzi.getTracks(cells,imdim);
    fprintf(1,'\t\t')
    toc

    %% Create the annotated tifs
    %   important for manual validation and checking takes about 20 sec per well, wells take about 40 sec total
    [cells2] = ProcessCzi.AnalyzeCells(cells2);
    % SC has uncommented this
    if DrawPlot
         tic
         ProcessCzi.DrawTracks(squeeze(im),cells2,wellTifSavePath,v)
         %ProcessCzi.DrawTracks2(squeeze(im),cells2,wellTifSavePath2);
         fprintf(1,'\t\t')
         toc
    end
    %% Save (x, y) coords in csv
    tic
    %ProcessCzi.ExportTrackStats(cells2,imdim,wellCsvSavePath)
    ProcessCzi.ExportTrackStatsSmall(cells2,imdim,wellCsvSavePath)
    fprintf(1,'\t\t')
    toc
    %         catch e
    %             fprintf(2,wellName)
    %             fprintf(2," exception: " + getReport(e)+"\n")
    %             s = {{[experPath,':']}, join({wellName,'failed'}), join({'exception:', getReport(e)})};
    %             % sendmail('brown.amelia@gmail.com','Update from matlab: RunProcessCzi',[s{:}])
    %             fprintf(2, "should have email\n")
    %             fprintf(2, datestr(datetime('now')) + "\n")
    %         end

end
%sendmail('brown.amelia@gmail.com','Matlab Finished: RunProcessCzi')

try
    %[exper,condits] = ReadCsvAsCondits(experPath);
    SummarizeData.RunSummarizeData(experPath)
catch

end

if ~(CONST.DRAW_MODE)
    close(f)
else
    close(v)
    close(v2)
end

end
