% % onset filter (roughly detect the note onsets)
% So = zeros(sizeS(1), sizeS(2));
% for i = 1:1:sizeS(1)
%     for j = 2:1:sizeS(2)
%         hi = Sg(i,j) - Sg(i,j-1);
%         if hi > 0.1
%             So(i,j) = hi;
%         end
%     end
% end
% figure;
% image(k,p,sfactor*So);
% set(gca,'YDir','normal');
% title('onset matrix');
% 
% % harmonic change filter (detect harmonic change boundaries)
% Sh = zeros(sizeS(1),sizeS(2));
% Shv = zeros(sizeS(1),sizeS(2)); % harmonic change matrix (one chord per col)
% Shc = zeros(1,sizeS(2));
% bassbound = 30;
% ht = 0.2;
% whs = 1;
% whe = 1;
% shidx = 1;
% for j = 1:1:sizeS(2)
%     for i = 1:1:bassbound
%         if So(i,j) > ht && j > 1 && j - whs > 10
%             % take the mean over the harmonic window in terms of row
%             whe = j-1;
%             wh = whs:whe;
%             for ii = 1:1:sizeS(1)
%                 gesiiwh = mean(Sg(ii,wh));
%                 if gesiiwh > 0.10
%                     Sh(ii,wh) = ones(1,length(wh))*gesiiwh;
%                 end
%             end
%             % normalize the content within harmonic window in terms of col
%             for jj = whs:1:whe
%                 tmp = Sh(:,jj);
%                 tmp = tmp / max(tmp);
%                 tmp(tmp < 0.2) = 0;
%                 Sh(:,jj) = tmp;
%             end
%             % fill the harmonic change vector
%             Shv(:,shidx) = Sh(:,whs);
%             Shc(shidx) = whs;
%             shidx = shidx + 1;
%             whs = j;
%             break;
%         end
%     end
% end
% nchords = shidx - 1;
% Shc(shidx) = sizeS(2);
% Shv = Shv(:,(1:nchords));
% Shc = Shc(1:nchords);
% figure;
% image(k,p,sfactor*Sh);
% set(gca,'YDir','normal');
% title('harmonic bounded salience matrix');
% figure;
% image(k(1:nchords),p,sfactor*Shv);
% set(gca,'YDir','normal');
% title('harmonic change matrix');