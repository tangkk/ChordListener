function plotchromagram(input, plotrange)

    sfactor = 100;
    figure;
    image(plotrange,1:length(input(:,1)),sfactor*input(:,plotrange));
    set(gca,'YDir','normal');
    title('input');
    xlabel('slice');
    ylabel('value');