function Fr = PlotOverCondition(Data,well,featureNames,treatStr,treatIdx,Groups,pageNum)
  %% Each Feature
  
  
  Pidx = 1;
  
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
            D2 = Data{i,2};
            
            D1 = D1(bCL,:);
            D2 = D2(bCL,:);
            
            Y2 = nanmean(D1,2);
            %         Slope = D2(:,1);
            %         Start = D2(:,2);
            %
            %         Y2 = Start + Slope*50;
            
            yrange = prctile(Y2,[0 99.5],1);
            yrange(2) = max(yrange(1)+0.1 , yrange(2));
            
            %% EACH TREATMENT
            for i2 = 1:length(treatStr)
                
                gidx =  W >= treatIdx(i2) &  W < treatIdx(i2+1);
                
                Slope1 = D2(gidx,1);
                Start1 = D2(gidx,2);
                Y22 = Y2(gidx);
                
                featname = featureNames{i,1};
                featUnit = featureNames{i,2};
                %% Start vs Growth
                subplot(1,3,1)
                
                title([titlestr, ': ', featname,' Increase over ' ,featname,' Start'])
                scatter(Start1,Slope1, 10, lineclrs(i2,:),'filled')
                xlabel([featname,' Start'])
                ylabel([featname,' Growth'])
                grid on
                hold on
                
                legend(treatStr,'location','southoutside')
                
                %% Over Concentration
                subplot(2,3,[2,3])
                title([featname,' Mean over Drug Concentration'])
                
                jitter = normrnd(0,.1, [1,nnz(gidx)] );
                
                scatter(Conc(gidx)+jitter, Y22 ,20,'filled')
                
                xlabel('Drug Concentration Step')
                ylabel([featname,' Mean ',featUnit])
                xlim([-1,13])
                ylim(yrange)
                grid on
                hold on
                
                %%
                %%subplot(2,3,4)
                
                %%
                subplot(2,3,[5,6])
                title([featname,' Mean over Drug Concentration'])
                
                for i3 = 0:max(Conc)
                    
                    X1 = i3;
                    X1 = X1 +  ((i2)/(length(treatStr)+1) -0.5)   /2;
                    
                    cidx = Conc(gidx)==i3;
                    if nnz(cidx)==0; continue; end
                    
                    Ysub = Y22(cidx);
                    Ymed = nanmedian(Ysub);
                    Yqrts = prctile(Ysub,[25 75],1);
                    %                 Ymean = mean(Ysub);
                    %                 Ysdev = std(Ysub);
                    
                    %% Medline Marker
                    X2 = X1 + [-0.1,0.1];
                    Y3 = Ymed*[1,1];
                    plot(X2,Y3, 'color', lineclrs(i2,:),'linewidth',3)
                    hold on
                    %% Vertical Lines
                    X2 = X1 + [0 0];
                    Y4 = Yqrts;
                    plot(X2,Y4,'color', lineclrs(i2,:),'linewidth',1)
                    
                end
                
                xlabel('Drug Concentration Step')
                ylabel([featname,' Mean'])
                xlim([-1,13])
                ylim(yrange)
                
                grid on
                hold on
                
            end
            
                Fr{Pidx} = getframe(gcf);
                Pidx = Pidx+1;
                pageNum = pageNum+1;

        end
    end

end 