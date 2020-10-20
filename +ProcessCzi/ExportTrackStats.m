%% ExportTrackStats  -  beta
%   from Epl
%   I've made a fair amount of changes
%


function ExportTrackStats(Cells,Dims,filePath)
    fprintf(1,'\tMaking csv files\n')

    %%
%     [P,F,~] = fileparts(Name);
%     P2 = fullfile(P,'Xls');
%     if ~exist(P2,'dir')
%         mkdir(P2)
%     end
%     XlsName = fullfile(P2,[F,'.csv']);

    %% Put Segs Together
    AllSegs = vertcat(Cells{:});
    AllTracks = [AllSegs.Tid];

    Tracks = unique(AllTracks);
    counts = hist(AllTracks,Tracks);
    Tracks = Tracks(counts>20);

    %%
    colCount = 5;
    %colCount = 4;

    Mat = cell(length(Cells)+colCount,colCount*max(length(Tracks),1));

    Mat(2,1:colCount:end) = {'Velocity(pixels/frame)'};
    Mat(2,2:colCount:end) = {'Tortuosity(pixels/pixels)'};
    Mat(2,3:colCount:end) = {'Mean Area(pixels.^2)'};
    Mat(2,4:colCount:end) = {'Mean Covarience'};

    Mat(4,1:colCount:end) = {'x(Pixel Position)'};
    Mat(4,2:colCount:end) = {'y(Pixel Position)'};
    Mat(4,3:colCount:end) = {'Area(pixels.^2)'};
    Mat(4,4:colCount:end) = {'Covariance'};
    Mat(4,5:colCount:end) = {'Pixels/Frame'};

    % imd = CONSTANTS.imageData;
    % ResolutionXYZ = imd.PixelPhysicalSize;
    % ResolutionTime = mean(diff(imd.TimeStampDelta)) * 60 * 60;

    for i = 1:length(Tracks)
        try
            Mat(1,(i-1)*colCount+1) = num2cell(Tracks(i));

            Tid = Tracks(i);

            cells = [AllSegs([AllSegs.Tid] == Tid)];

            %% Collect Size Stats
            Areas = [cells.Area]; %* ResolutionXYZ(1)*ResolutionXYZ(2);

            MeanSize = mean(Areas);

            %% Get Covariance
            HullsCR = {cells.PixelList};
            Hullscov = cellfun(@cov, HullsCR,'UniformOutput', false);
            try
                HullCov = cellfun(@(x) min(x(1,1),x(2,2)) / max(x(1,1),x(2,2)),Hullscov)';
            catch
                HullCov = nan(size(Hullscov,2),1);
            end
            %% Collect Movement Stats
            HullPos = vertcat(cells.Centroid);
            Times = [cells.time];
            [Times, I] = sort(Times);
            HullPos = HullPos(I,:);

            DeltaPos = [0,0; diff(HullPos,1,1)];
            Veloc = sqrt(DeltaPos(:,1).^2+DeltaPos(:,2).^2);
            MeanVoc = mean(Veloc);

            %Convert velocity to physical space and real time
            % MeanVoc = MeanVoc* ResolutionXYZ(1) / ResolutionTime;

            TotalDisp = sum(sqrt(DeltaPos(:,1).^2 + DeltaPos(:,2).^2));
            LengthTraveled = norm(HullPos(1,:) - HullPos(end,:));
            Tort =  TotalDisp./LengthTraveled;

            %% Write to Cell Mat


            Mat(3,(i-1)*colCount+1) = num2cell(MeanVoc);
            Mat(3,(i-1)*colCount+2) = num2cell(Tort);
            Mat(3,(i-1)*colCount+3) = num2cell(MeanSize);
            Mat(3,(i-1)*colCount+4) = num2cell(mean(HullCov));

            Mat(Times+4,((i-1)*colCount+1):((i-1)*colCount+2)) = num2cell(HullPos(:,[1,2]));
            Mat(Times+4,(i-1)*colCount+3) = num2cell(Areas);
            Mat(Times+4,(i-1)*colCount+4) = num2cell(HullCov);
            Mat(Times+4,(i-1)*colCount+5) = num2cell(Veloc);

        catch e

            % added because sometimes index exceeds matrix dimensions error
            % don't know what's causing it
            % it appears that for the same well it happens sometimes but not
            % always with what seems to be the same conditions
    %         disp("cell " + i + ", track " + Tracks(i))
    %         disp("exception: " + e.identifier)
    %         disp("except.message: " + e.message)
            fprintf(2,"cell " + i + ", track " + Tracks(i)+"\n")
            fprintf(2,"exception: " + getReport(e)+"\n")
            %keyboard
            continue

        end

    end
    T = cell2table(Mat);
    writetable(T,filePath,'WriteVariableNames',0);
    %disp(Mat)
    %xlswrite(XlsName,Mat)
end