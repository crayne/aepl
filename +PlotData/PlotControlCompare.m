function Fr = PlotControlCompare(Data,well,featureNames,treatStr,treatIdx,Groups,pageNum)


Pidx = 1;

%% Over Concentration with Significance
lineclrs = lines(6);
%% Summary Stat Figures
for g = 1:length(Groups)
    
    Grp = Groups{g};
    titlestr = [Grp, ' Stats'];
    bCL = strcmp(well(:,4),Grp);
    
    well1 = well(bCL,:);
    W = [well1{:,2}];
    W1 = [well1{:,3}];
    
    AllWells = {well1{:,5}};
    WellsUni = unique(AllWells);
    
    Conc = 14 - W1;
    Conc(W>= 17) = 0;
    
    %% Each Feature
    for i = 1:size(Data,1)
        
        %figure('position',[100,100,1800,800]);
        
        clf();
        AX1 = subplot(1,2,1);
        hold(AX1,'on')
        
        %%
        D1 = Data{i,1};
        D1 = D1(bCL,:);
        
        Y2 = nanmean(D1,2);
        
        yrange = prctile(Y2,[0 99.5],1);
        yrange(2) = max(yrange(1)+0.1 , yrange(2));
        
        %% EACH TREATMENT
        for i2 = length(treatStr) %% Control
            
            gidx =  W >= treatIdx(i2) &  W < treatIdx(i2+1);
            
            Y22 = Y2(gidx);
            
            featname = featureNames{i,1};
            featUnit = featureNames{i,2};
            
            Control = Y2(Conc(:)==0);
            
            
            %% Plot significance
            Xs = 1:length(WellsUni);
            Ys = nan(length(WellsUni),1);
            Ys2 = nan(length(WellsUni),1);

            WellsUni = unique(AllWells);
            WellNames = AllWells(gidx);

            for i3 = 1:length(WellsUni) %% For Each Well
                
                cidx = strcmp(WellNames,WellsUni{i3});
                
                Ysub = Y22(cidx);
                Ymed = nanmedian(Ysub);
                Yqrts = prctile(Ysub,[25 75],1);
                
                jitter = normrnd(0,.2, [1,nnz(cidx)] );
                scatter(i3+jitter, Ysub ,10,lineclrs(i3,:),'filled')   
                
                %% Medline Marker
                if i3==0
                    X2 = 0 + [-0.3,0.3];
                    X3 = [0 0];
                else
                    X2 = Xs(i3) + [-0.3,0.3];
                    X3 = Xs(i3) + [0 0];
                end
                Y3 = Ymed*[1,1];
                
                %%
                plot(AX1,X2,Y3, 'color', lineclrs(i3,:),'linewidth',5)
                hold on                
             
                
                %% Vertical Lines
                Yqrts = prctile(Ysub,[10 90],1);
                plot(AX1,X3,Yqrts,'color', lineclrs(i3,:),'linewidth',1)
                
                Yqrts = prctile(Ysub,[25 75],1);
                plot(AX1,X3,Yqrts,'color', lineclrs(i3,:),'linewidth',3)
                
                %% Vertical Lines
                Yqrts = prctile(Ysub,[45 55],1);
                plot(AX1,X3,Yqrts,'color', lineclrs(i3,:),'linewidth',8)

            end
        end
        
        title(AX1,[titlestr, ' : ' featname,' Control Values by Well'])
        xlabel(AX1,'Well ID Date')
        ylabel(AX1,[featname,' Control Values ', featUnit])
        
        xticks(AX1,1:length(WellsUni))
        xticklabels(AX1, WellsUni)
        xlim(AX1,[0,length(WellsUni)+1])
        ylim(AX1,yrange)
        
        grid(AX1,'on')
        AX1.Position = [0.05 0.1 0.9,0.85];
        
        Fr{Pidx} = getframe(gcf);
        Pidx = Pidx+1;
        pageNum = pageNum+1;
        
    end
end





end