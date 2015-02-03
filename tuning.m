% % tuning algorithm (this algorithm seems to be unreliable)
% Sbar = (sum(Spre,2))/sizeSpre(2); % first sum over all time frames
% fftSbar = fft(Sbar, 2^nextpow2(sizeSpre(1)));
% lenfftSbar = length(fftSbar);
% tuningpoint = lenfftSbar/6; % normalized frequency pi/3
% phi = angle(fftSbar(floor(tuningpoint)))*(tuningpoint - floor(tuningpoint))...
%     + angle(fftSbar(ceil(tuningpoint)))*(ceil(tuningpoint) - tuningpoint); % take the angle at pi/3
% d = - phi - 2*pi/3; % wrap d in [-pi,pi)
% if d < -pi
%     d = d + 2*pi;
% end
% if d > pi
%     d = d - 2*pi;
% end
% d = d / (2*pi);
% tau = 440*2^(d/12);
% if d < 0
%     for j = 1:1:sizeSpre(2)
%         col = Spre(:,j);
%         for i = 2:1:sizeSpre(1)
%             % do interpolation here
%             Spre(i,j) = -d*col(i-1) + (1+d)*col(i);
%         end
%     end
% else
%     for j = 1:1:sizeSpre(2)
%         col = Spre(:,j);
%         for i = 1:1:sizeSpre(1) - 1
%             % do interpolation here
%             Spre(i,j) = (1-d)*col(i) + d*col(i+1);
%         end
%     end
% end