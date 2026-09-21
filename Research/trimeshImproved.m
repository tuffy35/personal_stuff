function [output] = trimeshImproved(ele,x,y,fig)
    
    figure(fig)
    for i = 1:length(ele)
        ia = ele(i,1);
        ib = ele(i,2);
        ic = ele(i,3);

        xa = x(ia,1);
        ya = y(ia,1);

        xb = x(ib,1);
        yb = y(ib,1);

        xc = x(ic,1);
        yc = y(ic,1);

        xtemp = [xa xb xc xa];
        ytemp = [ya yb yc ya];
        
        hold on
        plot(xtemp,ytemp,'k-')
    end
    axis square
    xlim([0 1.01])
    output = 0;
end