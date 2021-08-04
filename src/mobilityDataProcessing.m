function mobilityDataProcessing()
%MOBILITYDATAPROCESSING process the NJIT 6DOF VR Navigation Dataset.
%
% References:
% M. Khan and J. Chakareski. NJIT 6DOF VR Navigation Dataset. 2020.
% https://www.jakov.org
%
% J. Chakareski and M. Khan, "Wifi-VLC dual connectivity streaming system
% for 6DOF multi-user virtual reality", in Proc. ACM Int'l Workshop on
% Network and Operating Systems Support for Digital Audio and Video
% (NOSSDAV), Istanbul, Turkey, Sep. 2021.
%
% NIST-developed software is provided by NIST as a public service. You may
% use, copy and distribute copies of the software in any medium, provided
% that you keep intact this entire notice. You may improve,modify and
% create derivative works of the software or any portion of the software,
% and you may copy and distribute such modifications or works. Modified
% works should carry a notice stating that you changed the software and
% should note the date and nature of any such change. Please explicitly
% acknowledge the National Institute of Standards and Technology as the
% source of the software. NIST-developed software is expressly provided
% "AS IS." NIST MAKES NO WARRANTY OF ANY KIND, EXPRESS, IMPLIED, IN FACT OR
% ARISING BY OPERATION OF LAW, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
% WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
% NON-INFRINGEMENT AND DATA ACCURACY. NIST NEITHER REPRESENTS NOR WARRANTS
% THAT THE OPERATION OF THE SOFTWARE WILL BE UNINTERRUPTED OR ERROR-FREE,
% OR THAT ANY DEFECTS WILL BE CORRECTED. NIST DOES NOT WARRANT OR MAKE ANY
% REPRESENTATIONS REGARDING THE USE OF THE SOFTWARE OR THE RESULTS THEREOF,
% INCLUDING BUT NOT LIMITED TO THE CORRECTNESS, ACCURACY, RELIABILITY,
% OR USEFULNESS OF THE SOFTWARE.
%
% You are solely responsible for determining the appropriateness of using
% and distributing the software and you assume all risks associated with
% its use,including but not limited to the risks and costs of program
% errors, compliance with applicable laws, damage to or loss of data,
% programs or equipment, and the unavailability or interruption of
% operation. This software is not intended to be used in any situation
% where a failure could cause risk of injury or damage to property.
% The software developed by NIST employees is not subject to copyright
% protection within the United States.
%
% 2020-2021 NIST/CTL (steve.blandino@nist.gov)



ts = 4e-3;
fftSize = 512;
measNum = 12;

yawId = 1; %1: YAW
rollId = 2; %2: ROLL
pitchId = 3;%3: PITCH

for i =1:measNum
    structLoad = load(sprintf('..\\Traces_6DOF_NJIT\\node%dmobility.mat',i));
    raw_data.(sprintf('MMs%d',i)) = structLoad.MMs1;
end

minTraceLen = min(structfun(@(x) size(x,1), raw_data));
nBlocks = fix(minTraceLen/fftSize);
traceLen = fftSize*nBlocks;
fs = 1/ts;
freqAx = linspace(-fs/2, fs/2-1/(2*fs), fftSize);

% Remove position information
rotationStruct = structfun(@removePositon, raw_data,'UniformOutput',false);
% Compute angular velocity
wStruct = structfun(@computeAngularVelocity, raw_data,'UniformOutput',false);

numFields  = length(fieldnames(rotationStruct));
traceNum = numFields;
rotationCell = cell(numFields,1);
wCell = cell(numFields,1);

for i = 1:numFields
    rotationCell{i} = rotationStruct.(sprintf('MMs%d', i))(1:traceLen,[1 3 2]);
    wCell{i} = wStruct.(sprintf('MMs%d', i));
end
rotationMat = cell2mat(rotationCell);
rotationMatRes = reshape(rotationMat, fftSize, nBlocks*traceNum,3);

% Create blocks for Welch's method
for i = 1: nBlocks*traceNum
    rotationStructChunks.(sprintf('MMs%d', i))= squeeze(rotationMatRes(:,i,:));
end

%% PLOT RAW DATA
figure(1)
subplot(3,1,1)
structfun(@(x) plotStruct(x, ts*(1:length(x)), yawId, 'Time (s)', ...
    'Orientation (deg)', 'YAW'), rotationStruct)
axis([0 120 -180 180])

subplot(3,1,2)
structfun(@(x) plotStruct(x, ts*(1:length(x)), pitchId, 'Time (s)', ...
    'Orientation (deg)', 'PITCH'), rotationStruct)
axis([0 120 -180 180])

subplot(3,1,3)
structfun(@(x) plotStruct(x, ts*(1:length(x)), rollId, 'Time (s)', ...
    'Orientation  (deg)', 'ROLL'), rotationStruct)
axis([0 120 -180 180])

figure(2)
plot(ts*(1:traceLen) , rotationMat(1:traceLen,yawId), 'g.', 'MarkerSize', 12)
hold on
plot(ts*(1:traceLen) , rotationMat(1:traceLen,pitchId), 'r.', 'MarkerSize', 12)
plot(ts*(1:traceLen) , rotationMat(1:traceLen,rollId), 'b.', 'MarkerSize', 12)
axis([0 120 -180 180])
ax = gca;
ax.YTick = [-180 -90 0 90 180];
p1 = plot(inf,inf, 'g', 'LineWidth', 4);
p2 = plot(inf,inf, 'r', 'LineWidth', 4);
p3 = plot(inf,inf, 'b', 'LineWidth', 4);
legend([p1 p2 p3], 'Yaw', 'Pitch', 'Roll', 'Location',  'best')
xlabel('Time (s)')
ylabel('Orientation  (deg)')
grid on

%% PLOT RAW HIST
%  Pitch
figure(4)
h = histc(rotationMat(:,pitchId), -180:0.01:180); %#ok<*HISTC>
[f,~,~] = fit((-180:0.01:180).', h/sum(h), 'gauss1');
fprintf('PDF Pitch: mu:%f, sigma:%f\n', f.b1, f.c1/sqrt(2));
p= plot(f, (-180:0.01:180).', h/sum(h));
p(2).LineWidth =4;
axis([-40 40 0 17e-4])
xlabel('Pitch (deg)');
ylabel('$f_\beta$ ','interpreter','latex');
grid on

% ROLL
figure(5)
h = histc(rotationMat(:,rollId), -180:0.01:180);
f = fit((-180:0.01:180).', h/sum(h), 'gauss1');
fprintf('PDF Roll: mu:%f, sigma:%f\n', f.b1, f.c1/sqrt(2));

p= plot(f, (-180:0.01:180).', h/sum(h));
p(2).LineWidth =4;
axis([-40 40 0 34e-4])
xlabel('Roll (deg)');
ylabel('$f_\gamma$ ','interpreter','latex');
grid on
hold on
p1 = plot(inf,inf, 'b.', 'MarkerSize', 24);
p2 = plot(inf,inf, 'r', 'LineWidth', 4);
legend([p1 p2], 'Measurements', 'Fitted Curve')

%% Plot clustered hist YAW

[id, C] = kmeansCenter(rotationMat(:,yawId),[-135 135 45 -45]);
id(id==max(unique(id))) = 1;
binSize = 1;
axHist =  (-180:binSize:180-binSize).';

for yawClusterId = 1: max(unique(id))
    yawCluster=C(id==yawClusterId);
    if yawClusterId ==1
        h = fftshift(histc(yawCluster, axHist));
        f = fit(axHist+180, h/sum(h), 'gauss1');
        fprintf('PDF Yaw%d: mu:%f, sigma:%f\n',yawClusterId, f.b1, f.c1/sqrt(2));
        figure(3)
        plot(axHist, fftshift(h/sum(h)), 'b.')
        hold on
        plot(axHist, [f(axHist(end/2:end)+180);f(axHist(1:end/2-1)+180)], 'r', 'LineWidth', 4)
    else
        
        yawCluster=C(id==yawClusterId);
        h = histc(yawCluster, axHist);
        f = fit(axHist, h/sum(h), 'gauss1');
        fprintf('PDF Yaw%d: mu:%f, sigma:%f\n',yawClusterId, f.b1, f.c1/sqrt(2));
        p= plot(f, axHist, h/sum(h));
        p(2).LineWidth =4;
        
    end
end
xlabel('Yaw (deg)');
ylabel('$f_\alpha$ ','interpreter','latex');
axis([-180 180 0 4e-2])
ax = gca;
ax.XTick = [-180 -90 0 90 180];
ax.YAxis.Exponent = -2;
grid on
legend off

%% Hist angular velocity
wStruct = cell2mat(wCell);
wYaw = reshape([wStruct.yaw]/ts, [],1);
wPitch = reshape([wStruct.pitch]/ts, [],1);
wRoll = reshape([wStruct.roll]/ts, [],1);

figure(7)
[h1, x] = hist(wYaw, fftSize,  'Normalization', 'pdf'); %#ok<*HIST>
p1= plot(x, h1/max(h1), 'Color', 'g', 'LineWidth', 4);
hold on
[h1, x] = hist(wPitch, fftSize,  'Normalization', 'probability');
p2 =  plot(x, h1/max(h1), 'Color', 'r', 'LineWidth', 4);
hold on
[h1, x] = hist(wRoll, fftSize,  'Normalization', 'probability');
p3 = plot(x, h1/max(h1), 'Color', 'b', 'LineWidth', 4);
legend([p1, p2, p3], 'Yaw', 'Pitch', 'Roll')
grid on
axis([-180 180 0 1])
xlabel('\omega (deg/s)')
set(gca, 'XTick', [-180 -135 -90 -45 0 45 90 135 180])

%% FFT
rotationStruct =rotationStructChunks;
rotationFFTStruct = structfun(@(x) fftStruct(x,fftSize), rotationStruct,'UniformOutput',false);

%% Plot PSD
psdYaw = structfun(@(x) plotStructSemiLog(x, freqAx, fs, yawId, 'Frequency (Hz)', 'PSD (dB/Hz)', 'YAW', 0),...
    rotationFFTStruct, 'UniformOutput' , false);
psdPitch = structfun(@(x) plotStructSemiLog(x, freqAx,  fs, pitchId, 'Frequency (Hz)', 'PSD (dB/Hz)', 'PITCH', 0),...
    rotationFFTStruct, 'UniformOutput' , false);
psdRoll = structfun(@(x) plotStructSemiLog(x, freqAx,  fs, rollId, 'Frequency (Hz)', 'PSD (dB/Hz)', 'ROLL', 0),...
    rotationFFTStruct, 'UniformOutput' , false);

%% Plot hist per subcarrier PSD and fit
% Plot Yaw PSD
psdBins= linspace(-134.5,24.5,50);
psdYawMat = cell2mat(struct2cell(psdYaw).');
for f = 1:fftSize
    yHistValue(:,f) = histcounts(psdYawMat(f,:).', 50,'BinLimits', [floor(psdBins(1)), ceil(psdBins(end))]);
end
yHistValue=yHistValue/max(yHistValue(:));
yHistValue(yHistValue == 0) =NaN;
[X,Y] = meshgrid(freqAx,psdBins);
figure(9)
h=pcolor(X,Y, yHistValue);
set(h, 'EdgeColor', 'none');
xlabel('$\nu$ (Hz)','interpreter','latex')
ylabel('$S_\alpha$ (dB/Hz)','interpreter','latex')
grid on

% Get median PSD
[~,i] = max(yHistValue);

toFit = psdBins(i);
toFit(1:fftSize/2) = [];
f = fit((1:150).', toFit(1:150).', 'exp2');
fprintf('PSD Yaw: a:%f, b:%f, c:%f, d:%f\n', f.a, f.b, f.c, f.d);

figure(10)
h=pcolor(X,Y, yHistValue);
set(h, 'EdgeColor', 'none');
xlabel('$\nu$ (Hz)','interpreter','latex')
ylabel('$S_\alpha$ (dB/Hz)','interpreter','latex')
set(h, 'EdgeColor', 'none');
set(h, 'FaceAlpha', 0.5);
hold on
plot(freqAx(257:256+143), f(1:143).','r', 'LineWidth', 4)
plot(freqAx(256:-1:256-142), f(1:143).', 'r', 'LineWidth', 4)
plot(freqAx([256+143 512]), -[64.785 64.785 ], 'r', 'LineWidth', 4)
plot(freqAx([1 256-142]), -[64.785 64.785 ], 'r', 'LineWidth', 4)
grid on

% Plot Pitch PSD
clear yHistValue
psdPitchMat = cell2mat(struct2cell(psdPitch).');
for f = 1:fftSize
    yHistValue(:,f) = histcounts(psdPitchMat(f,:).', 50,'BinLimits', [-135 25]);
end
yHistValue=yHistValue/max(yHistValue(:));
yHistValue(yHistValue == 0) =NaN;

[X,Y] = meshgrid(freqAx,psdBins);
figure(11)
h=pcolor(X,Y, yHistValue);
set(h, 'EdgeColor', 'none');
xlabel('$\nu$ (Hz)','interpreter','latex')
ylabel('$S_\beta$ (dB/Hz)','interpreter','latex')
grid on

[~,i] = max(yHistValue);
toFit = psdBins(i);
toFit(1:fftSize/2) = [];
f = fit((1:256).', toFit.', 'exp2');
fprintf('PSD Pitch: a:%f, b:%f, c:%f, d:%f\n', f.a, f.b, f.c, f.d);

figure(12)
h=pcolor(X,Y, yHistValue);
set(h, 'EdgeColor', 'none');
set(h, 'FaceAlpha', 0.5);
xlabel('$\nu$ (Hz)','interpreter','latex')
ylabel('$S_\beta$ (dB/Hz)','interpreter','latex')
grid on
hold on
plot(freqAx, [flip(f(1:256)); f(0:255)], 'r', 'LineWidth', 4)

% Plot Roll PSD
clear yHistValue
psdRollMat = cell2mat(struct2cell(psdRoll).');
for f = 1:fftSize
    yHistValue(:,f) = histcounts(psdRollMat(f,:).', 50,'BinLimits', [-135 25]);
end
yHistValue=yHistValue/max(yHistValue(:));
yHistValue(yHistValue == 0) =NaN;

[X,Y] = meshgrid(freqAx,psdBins);
figure(13)
h=pcolor(X,Y, yHistValue);
set(h, 'EdgeColor', 'none');
xlabel('$\nu$ (Hz)','interpreter','latex')
ylabel('$S_\gamma$ (dB/Hz)','interpreter','latex')
grid on

[~,i] = max(yHistValue);
toFit = psdBins(i);
toFit(1:fftSize/2) = [];
f = fit((1:256).', toFit.', 'exp2');
fprintf('PSD Roll: a:%f, b:%f, c:%f, d:%f\n', f.a, f.b, f.c, f.d)
figure(14)
h=pcolor(X,Y, yHistValue);
set(h, 'EdgeColor', 'none');
set(h, 'FaceAlpha', 0.5);
xlabel('$\nu$ (Hz)','interpreter','latex')
ylabel('$S_\gamma$ (dB/Hz)','interpreter','latex')
grid on
hold on
plot(freqAx, [flip(f(1:256)); f(0:255)], 'r', 'LineWidth',4)

end

function x =removePositon(x)
x(:,1:3) = [];
end

function y = fftStruct(x,fftSize)
y = fftshift(fft(x,fftSize,1))/sqrt(fftSize);
end

function yplot= plotStructSemiLog(y,x,Fs, rotIndex, xlabelValue,ylabelValue, titleValue, varargin)
if isempty(varargin)
    isPlot = 1;
else
    isPlot = varargin{1};
end

yplot = 20*log10(abs(y(:, rotIndex)./Fs));
if isPlot
    hold on
    plot(x,yplot)
end
% plot(x,angle(y(:, rotIndex)))
xlabel(xlabelValue)
ylabel(ylabelValue)
title(titleValue)
grid on
end

function plotStruct(y,x,rotIndex, xlabelValue,ylabelValue, titleValue)
hold on
plot(x,y(:, rotIndex), 'Color', [60 60 59]/255)
xlabel(xlabelValue)
ylabel(ylabelValue)
title(titleValue)
grid on
end

function w = computeAngularVelocity(x)
w.yaw = unwrap(x(2:end, 4)) - unwrap(x(1:end-1, 4));
w.pitch = x(2:end, 5) - x(1:end-1, 5) ;
w.roll = x(2:end, 6) - x(1:end-1, 6) ;
end

function  [id, C ] = kmeansCenter(x, center)

center = sort(center(:));
x = x(:);
lnCenter = length(center);
lnX = length(x);
id = zeros(lnX,1);
C = zeros(lnX,1);

for i = 1:lnCenter+1
    if i ==1
        idI = x<=center(i);
    elseif i== lnCenter+1
        idI = x>center(i-1);
    else
        idI = x<=center(i) & x>center(i-1);
    end
    id(idI) = i;
    C(idI) = x(idI);
end
end