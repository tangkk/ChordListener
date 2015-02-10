
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