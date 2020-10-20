function Fr = PlotDrugBreakdown(Data,well,featureNames,treatStr,treatIdx,Groups,pageNum)

Pidx = 1;

%% EACH GROUP
for g = 1:length(Groups)
    
    Grp = Groups{g};
    titlestr = [Grp, ' Stats'];
    bCL = strcmp(well(:,4),Grp);
    
    well1 = well(bCL,:);
    W = [well1{:,2}];
    W1 = [well1{:,3}];
    
    Conc = 14 - W1;
    Conc(W>= 17) = 0;
    lineclrs = jet(max(Conc)+1);
    
    %% EACH FEATURE
    
    for i = [3,4,5] %1:size(Data,1)
        
        %%
        D1 = Data{i,1};
        D1 = D1(bCL,:);
        
        %% EACH TREATMENT
        %             figure('position',[100,100,1800,800]);
        
        nTreatments = length(treatStr);
        nConc = max(Conc);
        idx1 = 1;
             
        %% Plot
        for i2 = 1:nTreatments
            
            clf();
            
            gidx =  W >= treatIdx(i2) &  W < treatIdx(i2+1);
            
            Y = D1(gidx,:)';
            ymax = max(Y(:))*1.1;
            ymax = max(ymax,0.1);
            featname = featureNames{i,1};
            featUnit = featureNames{i,2};
            Treat = treatStr{i2};
            
            idx1 = 1;
            
            AX1 = axes('Position',[0.45,0.97,0.1,0.1],'units','normalized');
            axis off;
            text(AX1,0,0,[Grp, ': ' ,Treat],'units','normalized','fontsize',20)
            
            for i3 = 0:(nConc)
                
                subplot(4,3,idx1)
                
                cidx = Conc(gidx)==i3;
                if nnz(cidx)==0; continue; end
                
                plot(Y(:,cidx));
                %plot(Y(:,cidx),  'color',[i3/nConc, 0, 0])
                hold on
                
                title(['Concentration Step ', num2str(i3),':',featname,' over Time'])
                xlabel('Time')
                ylabel([featname,' ',featUnit])
                ylim([-.1,ymax]);
                grid on
                hold on

                idx1 = idx1+1;
                
                AX5 = axes('Position',[0.95,0.05,0.1,0.1],'units','normalized');
                axis off;
                text(AX5,0,0,num2str(pageNum),'units','normalized','fontsize',20)
            end
            
            Fr{Pidx} = getframe(gcf);
            Pidx = Pidx+1;
            pageNum = pageNum+1;
        end
    end
end

end

