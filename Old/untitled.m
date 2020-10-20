    for i = 1:6   
        h = handles(i);
        set(h, 'PaperUnits','inches','PaperPosition',[0 0 11 8.5],'PaperOrientation','landscape')
        
        
        fname = [num2str(h.Number), '.pdf'];
        
        print(h,fname,'-dpdf')
        %close(h)
        
    end