function [seriesMetadata, varargout] = GetMetadata( bfReader, datasetExt )
%GETMETADATA Summary of this function goes here
%   Detailed explanation goes here

seriesMetadata = {};

if (~exist('datasetExt','var') || isempty(datasetExt))
    datasetExt = '';
end

orgMetadata = bfReader.getSeriesMetadata();
omeMetadata = bfReader.getMetadataStore();

onlyOneSeries = true;
if (bfReader.getSeriesCount()>1)
    prgs = Utils.CmdlnProgress(bfReader.getSeriesCount(),true);
    onlyOneSeries = false;
end

for series=0:bfReader.getSeriesCount()-1
    bfReader.setSeries(series);

    imageData = [];
    
    imageData.DatasetName = parseImageName(char(omeMetadata.getImageName(series)), datasetExt, bfReader.getSeriesCount());

    imageData.Dimensions = [bfReader.getSizeX(),bfReader.getSizeY(),bfReader.getSizeZ()];

    imageData.NumberOfChannels = double(omeMetadata.getChannelCount(series));
    imageData.NumberOfFrames = bfReader.getSizeT();

    xPixelPhysicalSize = safeGetValue(omeMetadata.getPixelsPhysicalSizeX(series));
    if xPixelPhysicalSize==0
        xPixelPhysicalSize = 1;
    end

    yPixelPhysicalSize = safeGetValue(omeMetadata.getPixelsPhysicalSizeY(series));
    if yPixelPhysicalSize==0
        yPixelPhysicalSize = 1;
    end

    zPixelPhysicalSize = safeGetValue(omeMetadata.getPixelsPhysicalSizeZ(series));
    if zPixelPhysicalSize==0
        zPixelPhysicalSize = 1;
    end

    imageData.PixelPhysicalSize = double([xPixelPhysicalSize, yPixelPhysicalSize, zPixelPhysicalSize]);

    if (strcmp(datasetExt,'.czi'))
        imageData.Position = [orgMetadata.get('Global Information|Image|S|Scene|Position|X #1'),...
                              orgMetadata.get('Global Information|Image|S|Scene|Position|Y #1'),...
                              orgMetadata.get('Global Information|Image|S|Scene|Position|Z #1')];
    elseif (omeMetadata.getPlaneCount(series)>0)
        imageData.Position = [safeGetValue(omeMetadata.getPlanePositionX(series,0)),...
                              safeGetValue(omeMetadata.getPlanePositionY(series,0)),...
                              safeGetValue(omeMetadata.getPlanePositionZ(series,0))];
    end

    imageData.ChannelNames = cell(imageData.NumberOfChannels,1);
    for c=1:imageData.NumberOfChannels
        colr = deblank(char(omeMetadata.getChannelName(series,c-1)));

        if (isempty(colr))
            colr = sprintf('Channel:%d',c);
        end

        imageData.ChannelNames{c} = colr;
    end

    imageData.StartCaptureDate = char(omeMetadata.getImageAcquisitionDate(series));
    ind = strfind(imageData.StartCaptureDate,'T');
    if (~isempty(ind))
        imageData.StartCaptureDate(ind) = ' ';
    end

    imageData.TimeStampDelta = 0;

    order = char(omeMetadata.getPixelsDimensionOrder(series));

    if (onlyOneSeries)
        prgs = Utils.CmdlnProgress(imageData.NumberOfFrames*imageData.NumberOfChannels*imageData.Dimensions(3),true);
        i = 1;
    end
    
    for t=1:imageData.NumberOfFrames
        for z=1:imageData.Dimensions(3)
            for c=1:imageData.NumberOfChannels
                ind = calcPlaneInd(order,z,c,t,imageData);
                try
                    delta = omeMetadata.getPlaneDeltaT(series,ind-1);
                catch er
                    delta = [];
                end
                if (~isempty(delta))
                    imageData.TimeStampDelta(z,c,t) = double(delta.value);
                end
                if (onlyOneSeries)
                    prgs.PrintProgress(i);
                    i = i+1;
                end
            end
        end
    end

    if (size(imageData.TimeStampDelta,1)~=imageData.Dimensions(3) ||...
            size(imageData.TimeStampDelta,2)~=imageData.NumberOfChannels || ...
            size(imageData.TimeStampDelta,3)~=imageData.NumberOfFrames)
        imageData = rmfield(imageData,'TimeStampDelta');
    end
    
    imageData.imageDir = fileparts(char(bfReader.getCurrentFile));

    seriesMetadata{series+1} = imageData;

    prgs.PrintProgress(series+1);
end

prgs.ClearProgress();

if (length(seriesMetadata)==1)
    seriesMetadata = seriesMetadata{1};
end

if (nargout>1)
    varargout{1} = omeMetadata;
end
if (nargout>2)
    varargout{2} = orgMetadata;
end
end

% Remove extension from image names while maintaining series information
function datasetName = parseImageName(imageName, datasetExt, numSeries)
    datasetName = imageName;

    extPattern = '(\.\w+?)';
    if ( ~isempty(datasetExt) )
        extPattern = ['(' regexptranslate('escape', datasetExt) ')'];
    end
    
    tokMatch = regexp(imageName,['(.+?)' extPattern '\s*(.*)'], 'tokens','once');
    if ( isempty(tokMatch) )
        return;
    end
    
    % Drop the series suffix if there's only one
    if ( numSeries == 1 )
        datasetName = tokMatch{1};
        return;
    end
    
    datasetName = [tokMatch{1} '_' parseSuffix(tokMatch{3})];
end

function seriesSuffix = parseSuffix(suffixString)
    seriesSuffix = suffixString;
    
    % Strip leading spaces and surrounding parens
    tokMatch = regexp(suffixString,'^\s*\(?(.+?)\)?\s*$','tokens','once');
    if ( isempty(tokMatch) )
        return;
    end
    
    seriesSuffix = tokMatch{1};
    
    abbrev = {'series','s'; 'position','p'};
    for i=1:size(abbrev,1)
        seriesSuffix = regexprep(seriesSuffix,[abbrev{i,1} '\s*' '(\d+)'], [abbrev{i,2} '$1']);
    end
end

function ind = calcPlaneInd(order,z,c,t,imageData)
switch order(3)
    case 'Z'
        ind = z-1;
        mul = imageData.Dimensions(3);
    case 'C'
        ind = c-1;
        mul = imageData.NumberOfChannels;
    case 'T'
        ind = t-1;
        mul = imageData.NumberOfFrames;
end

switch order(4)
    case 'Z'
        ind = ind + (z-1)*mul;
        mul = imageData.Dimensions(3)*mul;
    case 'C'
        ind = ind + (c-1)*mul;
        mul = imageData.NumberOfChannels*mul;
    case 'T'
        ind = ind + (t-1)*mul;
        mul = imageData.NumberOfFrames*mul;
end

switch order(5)
    case 'Z'
        ind = ind + (z-1)*mul;
    case 'C'
        ind = ind + (c-1)*mul;
    case 'T'
        ind = ind + (t-1)*mul;
end

ind = ind +1;
end

function val = safeGetValue(varIn)
if (isempty(varIn))
    val = 0;
    return
end

val = double(varIn.value);
end
