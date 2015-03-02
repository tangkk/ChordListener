function myLinePlot(x, y, xlab, ylab, xl, yl, marker, tit)

figure;
plot(x,y,marker);
xlim([1 xl]);
ylim([1 yl]);
xlabel(xlab);
ylabel(ylab);
title(tit);