% some old codes here that might be useful

% % calculate spectral centroid
% Sec = S(:,1:100);
% SX = sum(Sec,2);
% sc = round(sum(SX.*(1:length(SX))') / sum(SX));
% scw = sc/(length(SX));

% rmeanS = zeros(sizeSpre);
% rstdS = zeros(sizeSpre);
% rmeanC = zeros(sizeSpre);
% rstdC = zeros(sizeSpre);
% for j = 1:1:sizeSpre(2)
%     % do a running mean and std of the 216 bins within a sliding window of <18 bin
%     colS = Ss(:,j);
%     colC = Sc(:,j);
%     for i = 1:1:sizeSpre(1)
%         wmean = max(i-8,1):min(i+9,sizeSpre(1));
%         rmeanS(i,j) = mean(colS(wmean));
%         rmeanC(i,j) = mean(colC(wmean));
% %         rstdS(i,j) = std(colS(wmean)); % TODO: how can I include this?
% %         rstdC(i,j) = std(colC(wmean));
%     end
% end

% compute preliminary salience matrix
% Spre = Ss.*Sc;
% for i = 1:1:sizeSpre(1)
%     for j = 1:1:sizeSpre(2)
%         if Ss(i,j) < rmeanS(i,j) || Sc(i,j) < rmeanC(i,j)
%             Spre(i,j) = 0;
%         end
%     end
% end
% sfactor = 100;
% pp = 1:sizeNS(1);
% kk = 1:sizeNS(2);
% figure;
% image(kk,pp,sfactor*S);
% set(gca,'YDir','normal');
% title('preliminary salience matrix');
% tuning - have little effect for common commercial songs, add later
% my tuning consider only Sbar = (sum(Spre,2))/sizeSpre(2);
% and see if the location of peaks different from correct values

% S = zeros(sizeNS(1)/3, sizeNS(2));
% for i = 1:3:sizeNS(1)
%     for j = 1:1:sizeNS(2)
%         S((i+2)/3,j) = (Spre(i,j) + Spre(i+1,j) + Spre(i+2,j));
%     end
% end
% sizeS = size(S);
% for i = 1:1:sizeS(2)
%     S(:,i) = S(:,i) / max(S(:,i));
% end
% 
% sfactor = 100;
% p = 1:sizeS(1);
% k = 1:sizeS(2);
% figure;
% image(k,p,sfactor*S);
% set(gca,'YDir','normal');
% title('note salience matrix');


% initial gestalt filter (fill in signals that is automatically filled in by
% humans initially)
% wg = 10;
% Sg = zeros(sizeS(1), sizeS(2));
% gesc = 0;
% gesct = 1;
% gest = 0.1;
% for i = 1:1:sizeS(1)
%     for j = wg/2:wg/2:sizeS(2)-wg/2
%         gesval = mean(S(i,j-wg/2+1:j+wg/2));
%         if gesval > gest && gesc >= gesct
%             Sg(i,j-wg/2+1:j+wg/2) = gesval;
%             gesc = gesc + 1;
%         elseif gesval > gest
%             gesc = gesc + 1;
%         else
%             gesc = 0;
%         end
%     end
% end

% % feeback gestalt filter (fill in signals within harmonic change boundaries)
% Sg = zeros(sizeS(1), sizeS(2));
% gest = 0.1;
% for shcidx = 1:1:length(Shc) - 1
%     wg = Shc(shcidx):Shc(shcidx+1); % using harmonic change as cues
%     for i = 1:1:sizeS(1)
%         gesval = mean(S(i,wg));
%         if gesval > gest
%             Sg(i,wg) = gesval;
%         end
%     end
% end

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

% % harmonic change filter (detect harmonic change boundaries)
% Sh = zeros(sizeS(1),sizeS(2)); % averge out salience matrix via harmonic boundaries
% Shv = zeros(sizeS(1),sizeS(2)); % harmonic change matrix (one chord per col)
% Shc = zeros(1,sizeS(2));
% bassbound = 30;
% whs = 1;
% whe = 1;
% shidx = 1;
% prevbass = 0;
% trackidx = 1;
% firsttime = 1;
% Shc(shidx) = trackidx;
% for j = 1:1:sizeS(2)
%     for i = 1:1:bassbound
%         if Sg(i,j) > 0
%             if firsttime == 1
%                 prevbass = i;
%                 firsttime = 0;
%                 break;
%             end
%             curbass = i;
%             if curbass ~= prevbass || j == sizeS(2)
%                 lenHarm = j - trackidx;
%                 % average the harmonic bounded window
%                 for ii = 1:1:sizeS(1)
%                     Sh(ii,trackidx:j-1) = mean(Sg(ii,trackidx:j-1))*ones(1,lenHarm);
%                 end
%                 Shv(:,shidx) = Sh(:,trackidx);
%                 Shc(shidx) = trackidx;
%                 trackidx = j;
%                 shidx = shidx+1;
%             end
%             prevbass = curbass;
%             break;
%         end
%     end
% end
% nchords = shidx - 1;
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

% % if within a gestalt window ahead there's a non-zero bin, compensate the
% % blank in the middle
% wg = 20;
% Sg = zeros(sizeS(1), sizeS(2));
% for i = 1:1:sizeS(1)
%     trackidx = 1;
%     isblank = 0;
%     for j = 1:1:sizeS(2)
%         if S(i,j) > 0
%             % compensate the gestalt
%             if isblank == 1
%                 lenBlank = j - trackidx;
%                 if lenBlank <= wg
%                     Sg(i,trackidx:j-1) = S(i,trackidx)*ones(1,lenBlank);
%                 end
%             end
%             Sg(i,j) = S(i,j);
%             isblank = 0;
%             trackidx = j;
%         else
%             isblank = 1;
%         end
%     end
% end
% figure;
% image(k,p,sfactor*Sg);
% set(gca,'YDir','normal');
% title('note gestalt salience matrix - 1');
% % input from above, if a piece of salience is shorter than a gestalt window, ignore it
% wg = 10;
% for i = 1:1:sizeS(1)
%     trackidx = 1;
%     islight = 0;
%     for j = 1:1:sizeS(2)
%         if Sg(i,j) == 0
%             if islight == 1;
%                 lenLight = j - trackidx;
%                 if lenLight <= wg
%                     Sg(i,trackidx:j-1) = zeros(1,lenLight);
%                 end
%             end
%             trackidx = j;
%             islight = 0;
%         else
%             islight = 1;
%         end
%     end
% end
% figure;
% image(k,p,sfactor*Sg);
% set(gca,'YDir','normal');
% title('note gestalt salience matrix - 2');





% % input from above, if a piece of salience is shorter than a gestalt window, ignore it
% wg = 10;
% Sgneg = zeros(sizeS(1), sizeS(2));
% for i = 1:1:sizeS(1)
%     trackidx = 1;
%     islight = 0;
%     for j = 1:1:sizeS(2)
%         if S(i,j) == 0
%             if islight == 1;
%                 lenLight = j - trackidx;
%                 if lenLight <= wg
%                     Sgneg(i,trackidx:j-1) = zeros(1,lenLight);
%                 end
%             end
%             trackidx = j;
%             islight = 0;
%         else
%             Sgneg(i,j) = S(i,j);
%             islight = 1;
%         end
%     end
% end
% figure;
% image(k,p,sfactor*Sgneg);
% set(gca,'YDir','normal');
% title('note gestalt salience matrix - 1');
% 
% % if within a gestalt window ahead there's a non-zero bin, compensate the
% % blank in the middle
% wg = 20;
% Sgpos = zeros(sizeS(1), sizeS(2));
% for i = 1:1:sizeS(1)
%     trackidx = 1;
%     isblank = 0;
%     for j = 1:1:sizeS(2)
%         if Sgneg(i,j) > 0
%             % compensate the gestalt
%             if isblank == 1
%                 lenBlank = j - trackidx;
%                 if lenBlank <= wg
%                     Sgpos(i,trackidx:j-1) = Sgneg(i,trackidx)*ones(1,lenBlank);
%                 end
%             end
%             Sgpos(i,j) = Sgneg(i,j);
%             isblank = 0;
%             trackidx = j;
%         else
%             isblank = 1;
%         end
%     end
% end
% figure;
% image(k,p,sfactor*Sgpos);
% set(gca,'YDir','normal');
% title('note gestalt salience matrix - 2');
% 
% Sg = Sgpos; % gestalt salience matrix

%         % merge slash chords from nearest chords with same bass
%         if strcmp(ct,'maj/3')
%             if mod(pcb-4-1,12)+1 == cb && ~strcmp(pct,'maj/3') && ~strcmp(pct,'maj/5')
%                 chordogram{1,i-1} = cb;
%                 chordogram{2,i-1} = 'maj/3';
%                 chordogram{3,i-1} = num2bass(cb);
%                 yes = 1;
%                 continue;
%             end
%             if mod(ncb-4-1,12)+1 == cb && ~strcmp(nct,'maj/3') && ~strcmp(nct,'maj/5')
%                 chordogram{1,i+1} = cb;
%                 chordogram{2,i+1} = 'maj/3';
%                 chordogram{3,i+1} = num2bass(cb);
%                 yes = 1;
%                 continue;
%             end
%         end
%         if strcmp(ct,'maj/5')
%             if mod(pcb-7-1,12)+1 == cb && ~strcmp(pct,'maj/3') && ~strcmp(pct,'maj/5')
%                 chordogram{1,i-1} = cb;
%                 chordogram{2,i-1} = 'maj/5';
%                 chordogram{3,i-1} = num2bass(cb);
%                 yes = 1;
%                 continue;
%             end
%             if mod(ncb-7-1,12)+1 == cb && ~strcmp(nct,'maj/3') && ~strcmp(nct,'maj/5')
%                 chordogram{1,i+1} = cb;
%                 chordogram{2,i+1} = 'maj/5';
%                 chordogram{3,i+1} = num2bass(cb);
%                 yes = 1;
%                 continue;
%             end
%         end
%         % merge sus2 and sus4 into maj
%         if strcmp(ct,'sus2') || strcmp(ct,'sus4')
%             chordogram{2,i} = 'maj';
%             yes = 1;
%             continue;
%         end


% % *************** Test nnls chroma ************%
% nnlschroma = zeros(numtones, lenSlice);
% for i = 1:1:lenSlice
%     display(i);
%     nnlschroma(:,i) = lsqnonneg(Mc', X(:,i));
% end

% % *************** Test nnls chroma ************%
% % build the tone profiles for calculating note salience matrix
% % each sinusoidal at frequency 'ftone' is generated via sin(2*pi*n*f/fs)
% % try nnls method on Ss
% fmin = 27.5; % MIDI note 21
% fmax = 1661; % MIDI note 92
% fratio = 2^(1/36);
% numtonesn = 256;
% wln = 512;
% wn = hamming(wln);
% Msn = zeros(numtonesn, wln/2); % simple tone profiles
% Mcn = zeros(numtonesn, wln/2); % complex tone profiles
% % the true frequency of the tone is supposed to lie on bin notenum*3-1,
% % e.g. A4 is bin 49*3-1 = 146, C4 is bin 40*3-1 = 119 (note that notenum is
% % not midinum, note num is the index of the key on a piano with A0 = 1)
% for toneidx = 1:1:numtonesn
%     ftone = fmin*(fratio^(toneidx-2));
%     stone = sin(2*pi*(1:wln)*ftone/fs).*wn';
%     ctone = (0.9*sin(2*pi*(1:wln)*ftone/fs) + 0.9^2*sin(2*pi*(1:wln)*2*ftone/fs) + ...
%         0.9^3*sin(2*pi*(1:wln)*3*ftone/fs) + 0.9^4*sin(2*pi*(1:wln)*4*ftone/fs)).*wn';
%     
%     ffttone = abs(fft(stone));
%     ffttone = ffttone(1:wln/2);
%     ffttone = ffttone / norm(ffttone,2);
%     
%     fftctone = abs(fft(ctone));
%     fftctone = fftctone(1:wln/2);
%     fftctone = fftctone / norm(fftctone,2);
%     Msn(toneidx,:) = ffttone;
%     Mcn(toneidx,:) = fftctone;
% end
% nnlschroma = zeros(numtonesn, lenSlice);
% for i = 1:1:lenSlice
%     display(i);
%     nnlschroma(:,i) = lsqnonneg(Mcn', Ss(:,i));
% end
% sfactor = 1000;
% sizennls = size(nnlschroma);
% p = 1:sizennls(1);
% k = 1:sizennls(2);
% figure;
% image(k,p,sfactor*nnlschroma);
% set(gca,'YDir','normal');
% title('nnls chroma');
% ****** the computation is too slow.... ****** %

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

% % write as LRC files the true chord boundaries recognized manually
% target = 'haoting';
% path = strcat('../AudioSamples/',target,'/chordsections/');
% lrcpath = strcat('../AudioSamples/',target,'/');
% 
% files = sortFiles(path);
% bdrys = zeros(1,length(files) + 1);
% bdrys(1) = 0;
% for i = 1:1:length(files)
%     filepath = strcat(path,files(i).name);
%     af = audioinfo(filepath);
%     ns = af.TotalSamples;
%     fs = af.SampleRate;
%     bdrys(i+1) = bdrys(i) + ns;
% end
% bdrys = bdrys ./ fs;
% 
% lrc = strcat(lrcpath,strtok(files(1).name,'.'),'.gt.lrc');
% fw = fopen(lrc,'w');
% formatSpec1 = '%s';
% formatSpec2 = '%s\n';
% lenOut = length(bdrys) - 1;
% for i = 1:1:lenOut
%     sec = bdrys(i);
%     timestr = strcat(num2str(floor(sec/60)),':',num2str(mod(sec,60)));
%     s = strcat('[',timestr,']','N',num2str(i));
%     fprintf(fw, formatSpec2, s);
% end
% sec = bdrys(end);
% timestr = strcat(num2str(floor(sec/60)),':',num2str(mod(sec,60)));
% s = strcat('[',timestr,']');
% fprintf(fw, formatSpec1, s);
% fclose(fw);