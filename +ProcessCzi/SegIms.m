%%  SegIms  -  Epl
%       from Epl
%       little to no changes on my part
%
function [cellsOut,BW] = SegIms(im,DrawPlot,videoHandle)
fprintf(1,'\tSegmenting Images\n')
%%

cellsOut = [];
BW = [];
%% turn off usePar because of crash when usePar = 1;
usePar = 0;


cells = cell(1,size(im,5));

%% PreProcess Images Globally 
im1 = mat2gray(im);
imf = median(im1,5);
imf = imgaussfilt(imf,20);
im1 = im1 - imf;

imf = median(im1,5);
imD = mat2gray(im1 - imf);
im1 = im1 - (1-imD).*imf;
%%
%     im = permute(medfilt3(permute(im,[1 2 5 4 3]),[3,3,3]),[1 2 5 4 3]);
Fr = cell(size(im,5),1);

%% Segment Cells in Parallel
if usePar && ~DrawPlot
    
    parfor t = 1:size(im,5)        % for frames     [x,y,c,s,t]?
        %for t = randi(size(im,5))
        
        I1 = im1(:,:,1,1,t);
        CC = ProcessCzi.SegTexture_MSKCC(I1,t,DrawPlot);
        cells{t} = CC;
    end
    
else
    %% Segment cells Normally
    
    h = figure('position',[50 50 1000 1000]);
    
    BW = false(size(im1,1),size(im1,2),size(im1,5));
    for t = 1:size(im,5)        % for frames     [x,y,c,s,t]?
        %for t = randi(size(im,5))
        
        
        I1 = im1(:,:,1,1,t);
        [CC,imBW] = ProcessCzi.SegTexture_MSKCC(I1,t,DrawPlot);
        
        
        cells{t} = CC;
        
        BW(:,:,t) = imBW;

        drawnow
        F = getframe(gcf);
        F = F.cdata;
        Fr{t} = F;
        
    end
    close(h)
end

%% Prepare Cell Structure
cells = cells(~cellfun('isempty',cells));

cellsOut = vertcat(cells{:});

for i = 1:length(cellsOut)
    cellsOut(i).id = i;
    cellsOut(i).Tid = i;
end




%% Write Segmentations to Video 
if DrawPlot && exist('videoHandle','var') && ~isempty(videoHandle)
    for i = 1:length(Fr)
        writeVideo(videoHandle,Fr{i})
    end
end

end

