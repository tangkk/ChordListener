function visualizeChordProgression(outchordogram, outboundaries, tt)

if nargin == 2
    figure;
    hold on;
    Y = -10:0.1:10;
    for i = 1:1:length(outboundaries)
        X = outboundaries(1,i)*ones(size(Y));
        plot(X,Y);
    end
    hold off;
    div = (max(Y) - min(Y) - 1) / length(outchordogram);
    for i = 1:1:length(outchordogram)
        x = outboundaries(i);
        text(x,10 - i*div,outchordogram{i});
    end
    xlabel('slice');
    ylabel('chord');
    title('updated chordprogression vs. slice');
end

if nargin == 3
    figure;
    hold on;
    Y = -10:0.1:10;
    for i = 1:1:length(outboundaries)
        X = tt(outboundaries(1,i))*ones(size(Y));
        plot(X,Y);
    end
    hold off;
    div = (max(Y) - min(Y) - 1) / length(outchordogram);
    for i = 1:1:length(outchordogram)
        x = tt(outboundaries(i));
        text(x,10 - i*div,outchordogram{i});
    end
    xlabel('time');
    ylabel('chord');
    title('updated chordprogression vs. time');
end