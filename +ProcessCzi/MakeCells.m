function [CC] = MakeCells(imBW,t,I1,I2)

    [B] = bwboundaries(imBW,'noholes');
    
    CC = regionprops(imBW,I1,'centroid','area','PixelList','MeanIntensity');
    CC2 = regionprops(imBW,I2,'MeanIntensity');

    for ii = 1:length(CC)
        CC(ii).id = 0;
        CC(ii).Tid = 0;
        CC(ii).time = t;
        CC(ii).Bound = B{ii};
        CC(ii).MitoScore = CC(ii).MeanIntensity;
        CC(ii).Label = 0;
        CC(ii).wasSplit = round( CC2(ii).MeanIntensity);
    end

    CC = rmfield(CC,'MeanIntensity');
end 