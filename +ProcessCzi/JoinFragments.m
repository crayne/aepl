function SegsOut = JoinFragments(Segs,Edges)

Cells = [vertcat(Segs.id),vertcat(Segs.Tid),vertcat(Segs.time)];
TrackList = ([Segs.Tid])';
TrackData = zeros(max(TrackList),6);


for i = 1:max(TrackList)
    
    if ~any(TrackList==i); continue; end
    
    TrackData(i,1) = min([Segs(TrackList==i).time]);
    TrackData(i,2) = max([Segs(TrackList==i).time]);
    TrackData(i,3) = i;
end
TrackData(all(TrackData==0,2),:) = [];

maxTime = max(TrackData(:,2));

for t = 5:maxTime
    
    bEE = ismember(Cells(:,2), TrackData(TrackData(:,2)<maxTime-5,3));
    bLS = ismember(Cells(:,2), TrackData(TrackData(:,2)>5,3));
end 

EarlyEndCID = Cells(bEE,:);
LateStartCID = Cells(bLS,:);

TrackID = linkWithCosts2( EarlyEndCID,LateStartCID,TrackData,Edges);



%% Change the Track IDS

for i = 1:size(TrackID,1)
    
    TrA =  TrackID(i,1);
    TrBIdx =  find(TrackList == TrackID(i,2));
    
    for ii = 1:length(TrBIdx)
        Segs(TrBIdx(ii)).Tid = TrA;
    end
end 

SegsOut = Segs;



% % if UseCost
% %     
% %     %% Get All the corresponding Start Cells for the Tracks of Interest
% %     [EarlyEndCIds] = GetCell3(conn,EarlyEndTID,[t,t]);
% %     [LateStartCIds ] = GetCell3(conn,LateStartTID,TimeWindow);
% %     
% %     linkWithCosts2(conn, EarlyEndCIds,LateStartCIds,Tracks)
% %     
% % else
% %     
% %     %% Get All the corresponding Start Cells for the Tracks of Interest
% %     [EarlyEndS,ComsEE ] = GetCellALL(conn,EarlyEndTID,t);
% %     [LateStartS,COMSLS  ] = getCellsALL2(conn,LateStartTID,TimeWindow);
% %     
% %     %% Get The Center of Masses
% %     COMS = [ComsEE;COMSLS];
% %     linkWithOverlap(conn,EarlyEndS,LateStartS,COMS,range)
% % end



end

function linkWithCosts(conn, EarlyEndCID,LateStartCID,Tracks)
%% M by 4 Array with Props CID,TID,TIME,MaxRadius

ALLCIDS = [EarlyEndCID(:,1); LateStartCID(:,1)];
numString = sprintf('%.0f,' , ALLCIDS);
ALLCIDSn = numString(1:end-1);


cmd = ['SELECT cellID_src,cellID_dst,cost FROM tblDistCC ',...
    'WHERE cellID_src IN (' ALLCIDSn ')' ,...
    'AND cellID_dst IN (' ALLCIDSn ') group by cellID_src'];
CCDist = fetch(conn,cmd);
CCDist = cellfun(@double,CCDist,'uniformOutput',false);
CCDist = cell2mat(CCDist);

%% Sort By Costs
[~,ind] = sort(CCDist(:,3));
CCDist = CCDist(ind,:);

%% Retrack Connectable Objects
TakenTID = [];
for i = 1:size(CCDist,1)
    EarlierTrackID = [];
    LaterTrackID = [];
    
    SrcID = CCDist(i,1);
    DesID = CCDist(i,2);
    
    if  any(SrcID==EarlyEndCID(:,1),1)
        EarlierTrackID = EarlyEndCID(SrcID == EarlyEndCID(:,1),2);
        if any(DesID == LateStartCID(:,1),1)
            LaterTrackID = LateStartCID(DesID == LateStartCID(:,1),2);
        end
        
    elseif  any(SrcID==LateStartCID(:,1),1)
        EarlierTrackID =  LateStartCID(SrcID == LateStartCID(:,1),2);
        if any(DesID == EarlyEndCID(:,1),1)
            LaterTrackID = EarlyEndCID(DesID == EarlyEndCID(:,1),2);
        end
    end
    
    if isempty(EarlierTrackID) || isempty(LaterTrackID)
        continue
    end
    
    %Dont Track To Self
    if (EarlierTrackID) == LaterTrackID
        continue
    end
    
    %Dont retake Tracks
    if any(TakenTID==EarlierTrackID) || any(TakenTID==LaterTrackID)
        continue
    end
    
    %Dont add track If it doesnt Lengthen Track
    LaterEndTime = Tracks(Tracks(:,3) == LaterTrackID,2);
    EarlierEndTime = Tracks(Tracks(:,3) == EarlierTrackID,2);
    if EarlierEndTime >= LaterEndTime
        continue;
    end
    
    laterStartTime = Tracks(Tracks(:,3) == LaterTrackID,1);
    earlierStartTime = Tracks(Tracks(:,3) == EarlierTrackID,1);
    if laterStartTime <= earlierStartTime
        continue;
    end
    
    
    TakenTID = [TakenTID,EarlierTrackID,LaterTrackID];
    
    Tracking.ChangeTrack(conn, EarlierTrackID,LaterTrackID,EarlierEndTime)
end


end

function TrackIDs = linkWithCosts2( EarlyEndCID,LateStartCID,Tracks,Edges)
%% CID  M by 4 Array with Props CID,TID,TIME,MaxRadius

%% Tracks is Max,MIN,ID

%% Edges is SRC ID to DES ID

ALLCIDS = [EarlyEndCID(:,1); LateStartCID(:,1)];

[~,ind] = sort( EarlyEndCID(:,2) - EarlyEndCID(:,1));
EarlyEndCID = EarlyEndCID(ind,:);

[~,ia,~] = intersect(Tracks(:,3),LateStartCID(:,2));
LateTracks = Tracks(ia,:);

%% Sort by Length of Track
[~,ind] = sort( LateTracks(:,2) - LateTracks(:,1));
LateTracks = LateTracks(ind,:);

%% Retrack Connectable Objects
TrackIDs = [0,0];
TakenTID = [];

for i = 1:size(EarlyEndCID,1)
    EarlierTrackID = EarlyEndCID(i,2);
    
    MatchList =  Edges(Edges(:,1) == EarlyEndCID(i,1),2);
    MatchList =  [MatchList; Edges(Edges(:,2) == EarlyEndCID(i,1),1)];
    
    if any(EarlierTrackID == TrackIDs(:,1))
        continue;
    end 

    for j = 1:size(LateTracks,1)
        
        LaterTrackID = LateTracks(j,3);
        LaterCID = LateStartCID(LateStartCID(:,2) == LaterTrackID,1);
        
        if isempty(EarlierTrackID) || isempty(LaterTrackID)
            continue
        end
        
        %Dont Track To Self
        if (EarlierTrackID) == LaterTrackID
            continue
        end
        
        %Dont retake Tracks
        if  any(TakenTID==LaterTrackID)
            continue
        end
        
        %% Check If Cost Exists Between
        CCMatches = intersect(MatchList,LaterCID);
        if isempty(CCMatches)
            continue
        end
        
        %% Dont add track If it doesnt Lengthen Track
        LaterEndTime = LateTracks(j,2);
        EarlierEndTime = Tracks(Tracks(:,3) == EarlierTrackID,2);
        if EarlierEndTime >= LaterEndTime
            continue;
        end
        
        laterStartTime = LateTracks(j,1);
        earlierStartTime = Tracks(Tracks(:,3) == EarlierTrackID,1);
        if laterStartTime <= earlierStartTime
            continue;
        end
        
        %% Skip if its already been Tracked to 
        if any(LaterTrackID == TrackIDs(:,2))
            continue 
        end 
        
        TrackIDs = [TrackIDs; [EarlierTrackID,LaterTrackID]];
        %disp(TrackIDs)
        
        
        %Tracking.ChangeTrack(EarlierTrackID,LaterTrackID,EarlierEndTime)
        
        break;
        
    end
    
end

TrackIDs(1,:) = [];

end






