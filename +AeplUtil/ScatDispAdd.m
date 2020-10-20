function ScatDispAdd(data1,data2)


%     blu =   [0.00,0.45,0.75];
%     ora =   [0.85,0.33,0.10];
%     yel =   [0.93,0.69,0.13];
%     purp =  [0.49,0.18,0.56];
%     gree =  [0.47,0.67,0.19];
%     lblu =  [0.30,0.75,0.93];
%     red =   [0.64,0.08,0.18];
%     
%     myColors = [blu;ora;yel;gree;purp;lblu;red]
    
    if nargin == 1
        boop = size(data1);
        if boop(1) == 2
            x = data1(1,:);
            y = data1(2,:);
            count = boop(2);
        elseif boop(2) == 2
            x = data1(:,1);
            y = data1(:,2);
            count = boop(1);
        end
    elseif nargin == 2
        x = data1;
        y = data2;
        boop = size(data1);
        if boop(1) == 1
            count = boop(2);
        elseif boop(2) == 1
            count = boop(1);
        end
    end
    
    for i =1:count 
        y(i) = 1024-y(i);
    end
    
    %count = size(x)
    hold on
    cmap = AeplUtil.BuildCmap([1:count]);
    
    scatter(x,y,[],cmap,'filled')
    %get(c)
    hold off
    ylim([0,1024])
    xlim([0,1024])
    
    %hold on
%     while i < count+1
%     
%     %figure('visible','off');
%     
%         scatter(x[i],y,filled)
    
end