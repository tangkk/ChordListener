% plot various gram as an image with a scaling factor
function myImagePlot(img, x, y, xl, yl, tit, ytl, ytlab)
if nargin == 6
    sfactor = 100;
    figure;
    image(x,y,sfactor*img);
    set(gca,'YDir','normal');
    xlabel(xl);
    ylabel(yl);
    title(tit);
end
if nargin == 8
    sfactor = 100;
    figure;
    image(x,y,sfactor*img);
    set(gca,'YDir','normal');
    xlabel(xl);
    ylabel(yl);
    set(gca, 'YTick',ytl, 'YTickLabel', ytlab);
    title(tit);
end