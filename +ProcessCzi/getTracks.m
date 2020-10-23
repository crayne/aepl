%% GetTracks  - Epl
%   frpm epl
%   little to no changes
%


function [SegsOut,Edges] = getTracks(Segs,dims)

fprintf(1,'\tGetting Tracks from segmentation\n')

maxT = max([Segs.time]);

Edges = cell(1,maxT);

ids = vertcat(Segs.id);
times = vertcat(Segs.time);

TimeWindow = 10;
MaxDist = 150;

%% Get Edges 
%% Loop through frames
for t = 1:maxT-1
    %% Get all segments for time t (that is, for one frame)
    T0 = Segs(times==t);
    Trange1 = t+1;
    Trange2 = min(t+TimeWindow , length(Segs));
    %% For comparison, segments t+1 thru t+TimeWindow
    T1 = Segs(times>= Trange1 & times <= Trange2); 
    Edges{t} = MakeEdgeT(T0,T1,[ids,times],dims,MaxDist);
end

EdgesT = vertcat(Edges{:});
[~,idx] = sort(EdgesT(:,3),1);
SortT = EdgesT(idx,:);

SortT(SortT(:,1)==0 | SortT(:,2)==0,:) = [];

SegVect  = Segs;
TrackVect = [SegVect.id;SegVect.Tid]';
EdgeVect = NaN(length(Segs),1);

%% Assign Tracks
BTaken = false(length(TrackVect),1);
ATaken = false(length(TrackVect),1);

for i = 1:size(SortT,1)
    
    bDes = TrackVect(:,2)==SortT(i,2);
    bSrc = TrackVect(:,1)==SortT(i,1);
    
    if BTaken(SortT(i,2)) || ATaken(SortT(i,1))
        %if ismember(SortT(i,2),BTaken) || ismember(SortT(i,1),ATaken)
    else
        
        TrackVect(bDes,2) = TrackVect(bSrc,2);
        EdgeVect(bSrc) = SortT(i,3);
        %     ATaken = [ATaken,SortT(i,1)];
        %     BTaken = [BTaken,SortT(i,2)];
        
        ATaken(SortT(i,1)) = 1;
        BTaken(SortT(i,2)) = 1;
        
    end
end

%% 
% % for ii = 1:size(SortT,1)
% %
% %     id0 = SortT(ii,1);
% %     id1 = SortT(ii,2);
% %
% %     if id0 == 0; continue; end
% %
% %     if any(TakenA==id0) || any(TakenB==id1)
% %         continue
% %     end
% %
% %     if id1 < 0 || id0 < 0
% %         continue
% %     end
% %
% %     TakenA = [TakenA,id0];
% %     TakenB = [TakenB,id1];
% %
% %     TrkA = TrackVect(TrackVect(:,1)==id0,2);
% %     TrkB = TrackVect(TrackVect(:,1)==id1,2);
% %
% %     TrackVect(TrackVect(:,2)==TrkB,2) = TrkA;
% %
% % end



%% Prepare for output 
SegsOut = Segs;
for ii = 1:length(Segs)
    sid = Segs(ii).id;
    Tid = TrackVect(TrackVect(:,1)==sid,2);
    Eval = EdgeVect(TrackVect(:,1)==sid,1);
    
    SegsOut(ii).Edge = Eval;
    
    SegsOut(ii).Tid = Tid;
end

%% Mend broken tracks 
SegsOut = ProcessCzi.JoinFragments(SegsOut,EdgesT);



%% Plot 
if 0
    % plot
    f = figure;
    T = [SegsOut.Tid];
    for i = 1:max(T)
        bT = T==i;
        if nnz(bT)==0; continue; end
        pts = vertcat(SegsOut(bT).Centroid);
        plot(pts(:,1),pts(:,2))
        hold on
    end
    close(f)
end

end
%%

function CalcIOU = CalculateIOU(C1,C2)
Intersect = intersect(C1.PixelList, C2.PixelList, 'rows');
Union = union(C1.PixelList, C2.PixelList, 'rows');
CalcIOU = length(Intersect)/length(Union);
end


function EdgeList = MakeEdgeT(T0,T1,Tracks,dims,maxDist)
% Use Intersection Over Union tracking method
useIOU = 1;
AdjustCentroids = 0;
EdgeList = [0,0,0];
idx = 1;
for i = 1:length(T0)
    
    if any(EdgeList(:,1) == T0(i).id)
        continue
    end
    C1 = T0(i);
       
    StartTime = Tracks(Tracks(:,1)==C1.Tid,2);
    %if useIOU is true, go through once
    %if all the IOUs are 0, change useIOU to false
    %?Return useIOU because sort will go the other way
    AllZeroIOUs = 1;
    if useIOU == true
       for j = 1:length(T1)
           C2 = T1(j);
           IOU = CalculateIOU(C1,C2);
           if IOU ~= 0.0
               AllZeroIOUs = 0;
               break;
           end
       end
       if AllZeroIOUs == 1
           useIOU = false;
           AdjustCentroids = true;
           fprintf(1,'All IOUs here are zero -- using Centroids\n');
       end
    end
    for j = 1:length(T1)
        
        C2 = T1(j);
        %Use Intersection Over Union
        if useIOU == true
            IOU = CalculateIOU(C1,C2);
            if  IOU == 0.0
                continue;
            end
            % Add 1 to IOU so range is 1 < IOU < 2
            Edist = IOU + 1;
         else
            %Use Centroid 
            dist =  sqrt(sum((C1.Centroid-C2.Centroid).^2,2));
        
            if dist>maxDist
                continue;
            end 
            %I don't know what this code does, but Centroid doesn't work
            %without it
            dist = 2*dist/(C1.Area + C2.Area);       
            deltaT = abs(C1.time - C2.time)-1;       
            Edist = dist + deltaT;
            if (AdjustCentroids == true)
                Edist = -(Edist/maxDist);
            end                
 
        end
              
        %% Add Match to List
        EdgeList(idx,:) = [C1.id,C2.id,Edist];
        idx = idx+1;
    end
    
end

%% Calculate Escape Costs;

% for i = 1:length(T0)
%  C1 = T0(i);
%  dist1 = 2*C1.Centroid;
%  dist2 = 2*dims(1:2) - C1.Centroid;
%  dist = min([dist1,dist2]);
%  if dist > maxDist; continue; end
%  dist = dist/maxDist;
%  EdgeList(idx,:) = [C1.id,-1,dist];
%  idx = idx+1;
%
% end
%
% for i = 1:length(T1)
%  C2 = T1(i);
%  dist1 = 2*C2.Centroid;
%  dist2 = 2*dims(1:2) - C2.Centroid;
%  dist = min([dist1,dist2]);
%  if dist > maxDist; continue; end
%  dist = dist/maxDist;
%  EdgeList(idx,:) = [-1,C2.id,dist];
%  idx = idx+1;
% end

end

