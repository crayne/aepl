function Fr = PlotStatDiff(Data,well,featureNames,treatStr,treatIdx,Groups,pageNum)


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
    
    Conc = 14 - W1;
    Conc(W>= 17) = 0;
    
    %% Each Feature
    for i = 1:size(Data,1)
        
        %         figure('position',[100,100,1800,800]);
        clf();
        
        D1 = Data{i,1};
        D1 = D1(bCL,:);
        
        Y2 = nanmean(D1,2);
        
        yrange = prctile(Y2,[0 99.5],1);
        yrange(2) = max(yrange(1)+0.1 , yrange(2));
        
        AX1 = subplot(1,2,1);
        AX2 = subplot(1,2,2);
        hold(AX1,'on')
        %% EACH TREATMENT
        for i2 = 1:length(treatStr)
            
            gidx =  W >= treatIdx(i2) &  W < treatIdx(i2+1);
            
            Y22 = Y2(gidx);
            
            featname = featureNames{i,1};
            featUnit = featureNames{i,2};
            
            Control = Y2(Conc(:)==0);
            %% Plot significance
            Xs = nan(max(Conc),1);
            Ys = nan(max(Conc),1);
            Ys2 = nan(max(Conc),1);
            
            for i3 = 1:max(Conc)
                
                X1 = i3;
                X1 = X1 +  ((i2)/(length(treatStr)+1) -0.5)   *0.9;
                
                cidx = Conc(gidx)==i3;
                if nnz(cidx)==0; continue; end
                Ysub = Y22(cidx);
                
                [~,sigVal] = ttest2(Control,Ysub,'Vartype','unequal');
                Ys(i3) = sigVal;
                
                [p] = ranksum(Control,Ysub);
                Ys2(i3) = p;
                
                if 1 || (p < (10^-2) || sigVal <(10^-2))
                    %% Vertical Lines
                    %                 X2 = X1*ones(size(sigVal));
                    Xs(i3) = X1;
                end
            end
            
            %%
            
            plot(AX2,Xs,Ys2,'x-','color', lineclrs(i2,:),'linewidth',1)
            hold on
            plot(AX2,Xs,Ys,'o-','color', lineclrs(i2,:),'linewidth',1)
            
            hold on
            set(AX2, 'YScale', 'log')
            
            plot([2 10],[1,1]*10^-10,'k')
            plot([2 10],[1,1]*10^-15,'k')
            plot([2 10],[1,1]*10^-5,'k')
            grid on
            
            xlabel(AX2,'Drug Concentration Step')
            ylabel(AX2,'P value')
            title(AX2,'Statistical Difference from Control(T Test & Ranksum)')
            xlim(AX2,[-1,13])
            ylim(AX2,[10^-18,1])
            
            %legend(AX2,treatStr,'location','southoutside')
            AX2.Position = [0.65 0.2 0.3,0.75];
            %%
            
            for i3 = 0:max(Conc)
                
                cidx = Conc(gidx)==i3;
                if nnz(cidx)==0; continue; end
                
                Ysub = Y22(cidx);
                Ymed = nanmedian(Ysub);
                Yqrts = prctile(Ysub,[25 75],1);
                
                %% Medline Marker
                if i3==0
                    X2 = 0 + [-0.3,0.3];
                    X3 = [0 0];
                else
                    X2 = Xs(i3) + [-0.3,0.3];
                    X3 = Xs(i3) + [0 0];
                end
                Y3 = Ymed*[1,1];
                
                
                plot(AX1,X2,Y3, 'color', lineclrs(i2,:),'linewidth',5)
                hold on
                
                %% Vertical Lines
                Yqrts = prctile(Ysub,[10 90],1);
                plot(AX1,X3,Yqrts,'color', lineclrs(i2,:),'linewidth',1)
                
                Yqrts = prctile(Ysub,[25 75],1);
                plot(AX1,X3,Yqrts,'color', lineclrs(i2,:),'linewidth',3)
                
                %% Vertical Lines
                Yqrts = prctile(Ysub,[45 55],1);
                plot(AX1,X3,Yqrts,'color', lineclrs(i2,:),'linewidth',8)
            end
        end
        
        title(AX1,[titlestr, ' : ' featname,' Mean over Drug Concentration'])
        xlabel(AX1,'Drug Concentration Step')
        ylabel(AX1,[featname,' Mean ', featUnit])
        xlim(AX1,[-1,13])
        ylim(AX1,yrange)
        grid(AX1,'on')
        set(AX2, 'YDir','reverse')
        AX1.Position = [0.05 0.1 0.55,0.85];
        
        
        AX5 = axes('Position',[0.95,0.05,0.1,0.1],'units','normalized');
        axis off
        
        text(AX5,0,0,num2str(pageNum),'units','normalized','fontsize',20)
        
        Fr{Pidx} = getframe(gcf);
        Pidx = Pidx+1;
        pageNum = pageNum+1;
        
    end
end



end 