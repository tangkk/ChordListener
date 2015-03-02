% plot various gram as an image with a scaling factor
function myImagePlot(img, x, y, xl, yl, tit)

sfactor = 100;
figure;
image(x,y,sfactor*img);
set(gca,'YDir','normal');
xlabel(xl);
ylabel(yl);
title(tit);