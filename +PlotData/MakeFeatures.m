function Data = MakeFeatures(data)

%% Get Stats
DeathCountTime = zeros(length(data), max( data{1}.Time_FrameNum_));
MitoCountTime = zeros(length(data), max( data{1}.Time_FrameNum_));
PopulCountTime = zeros(length(data), max( data{1}.Time_FrameNum_));

for i = 1:length(data)
    d = data{i};
    
    L = d.Label;
    T = d.Time_FrameNum_;
    
    for t = 1:max(d.Time_FrameNum_)
        DeathCountTime(i,t) = nnz(L==2 & T==t);
        MitoCountTime(i,t) = nnz(L==1  & T==t);
        PopulCountTime(i,t) = nnz(L==0  & T==t);
        
    end
end

%% Tracks per stat

VelocityTime = zeros(length(data), max( data{i}.Time_FrameNum_));
PopulCountTime2 = zeros(length(data), max( data{i}.Time_FrameNum_));

% VeloCount = zeros(1, max( data{i}.Time_FrameNum_));
% VeloSum = zeros(1, max( data{i}.Time_FrameNum_));

maxTime = 97;

for i = 1:length(data)
    
    d = data{i};
    
    TIDs = [d.TrackID];
    CIDs = [d.CellID];
    Times = [d.Time_FrameNum_];
    X = [d.x_PixelPosition_];
    Y = [d.y_PixelPosition_];
    bSplit = [d.wasSplit];
    lState = [d.Label];
    
    VeloCount = zeros(1, maxTime);
    VeloSum = zeros(1, maxTime);
    
    for tid = 1:max(TIDs)
        
        bTr = TIDs == tid & bSplit==1 & lState==0;
        if nnz(bTr)<10; continue; end
        time = Times(bTr);
        locs = vertcat([X(bTr),Y(bTr)]);
        dist = sqrt(sum((diff(locs).^2)  ,2));
        delT = diff(Times(bTr));

        dist = dist ./ delT;
        time = time(1:end-1);
        
        VeloSum(time) = VeloSum(time) + dist';
        VeloCount(time) = VeloCount(time) +1;
        
    end
    

    PopCount = zeros(1, maxTime);
    for tid = 1:max(TIDs)
        bTr = TIDs == tid & lState==0;
        if nnz(bTr)<10; continue; end
        time = Times(bTr);
        PopCount(time) = PopCount(time)+1;
        
    end
        
    VelocityTime(i,:) = (VeloSum ./ VeloCount);
    PopulCountTime2(i,:) = PopCount;
end


%%
DeathCountTime = medfilt2(DeathCountTime,[1,5]);
MitoCountTime = medfilt2(MitoCountTime,[1,3]);
%PopulCountTime = medfilt2(PopulCountTime,[1,5]);
PopulCountTime2 = medfilt2(PopulCountTime2,[1,5]);
VelocityTime = medfilt2(VelocityTime,[1,5]);

VelocityTime = VelocityTime * (1/1.5504) * 10; 

Data = cell(4,2);
Data{1,1} = MitoCountTime;
Data{2,1} = DeathCountTime;

Data{3,1} = DeathCountTime./ (DeathCountTime + PopulCountTime2+.1);
% Data{3,1} = PopulCountTime;
Data{4,1} = PopulCountTime2;
Data{5,1} = VelocityTime;




%% Get Recursion Stats
for i = 1:size(Data,1)
    
    F = Data{i,1};    
    for ii = 1:size(F,1)
        
        
        d = F(ii,5:end-10);
        t = 1:length(d);
        t(isnan(d)) = [];
        d(isnan(d)) = [];
        
        p = polyfit(t,d,1);
        Data{i,2}(ii,1:2) = p;
        
    end
end



end 