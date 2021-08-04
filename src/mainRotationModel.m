%% Head Rotation Model for VR System Level Simulations.
% Steve Blandino, Tanguy Ropitault, Raied Caromi, Jacob Chakareski,
% Mahmudur Khan, and Nada Golmie. 2021.
%
% Script to generate VR head rotation angles (yaw, pitch and roll)

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

close all
%% Input Parameters
ts = 4e-3; % Sampling Time
yawCoh = 4096; % Yaw: samples generated with one distribution
sigmaN2 = 0.5;
isPlot = 1;


%% Dependent Params
fs = 1/ts; % Sampling Frequency
Ns = pow2(15); % Number of samples for ifft calculation
df = fs/Ns; % Frequency resolution
Nbb = Ns/2; % Equivalent baseband samples.
L = Ns/yawCoh; % Yaw: # of blocks (each block with different distribution)
freqAx = linspace(0, fs/2, Ns).'; % frequency index range

%% TMP

%% YAW PDF
yawMultiModalP = [0.23 0.4971 0.7551 1];
muYawMultiModal = [173, -92, -6, 87];
sigmaYawMultiModal = [28.36 25.20 24.78 26.61];

muYaw = zeros(L,1);
sigmaYaw = zeros(L,1);
yaw = zeros(yawCoh, L);

for i =1:L
    indexMultiModal = find((rand<yawMultiModalP) == 1, 1);
    muYaw(i) = muYawMultiModal(indexMultiModal);
    sigmaYaw(i) = sigmaYawMultiModal(indexMultiModal);
    yaw(:,i) = muYaw(i)+ sigmaYaw(i)*randn(yawCoh,1);
end

%% ROLL PDF
muRoll = 0.3186;
sigmaRoll = 4.048/sqrt(2);
roll = muRoll + sigmaRoll*randn(Ns,1);

%% PITCH PDF
muPitch = -3.083;
sigmaPitch = 15.37/sqrt(2);
pitch = muPitch + sigmaPitch*randn(Ns,1);

%% YAW
a = -49.15;
b = 0.001928;
c = 48.09;
d = -0.1768;
nu_c = 49.7147;

psd1side= biExpPsd(a,b,c,d,nu_c,freqAx);
yawModelOut = getRotation(yaw, psd1side, Ns, df, sigmaN2);

if isPlot
    figure
    subplot(3,1,1)
    plot((1:Ns)*4e-3, wrapTo180(yawModelOut), '*')
    ylim([-180 180])
end

%% ROLL
a = -35.65;
b = 0.001234;
c = 38.85;
d = -0.07319;
nu_c = fs/2;

psd1side= biExpPsd(a,b,c,d,nu_c,freqAx);
rollModelOut = getRotation(roll, psd1side, Ns, df);

if isPlot
    subplot(3,1,3)
    plot((1:Ns)*4e-3, wrapTo180(rollModelOut), '*')
end

%% PITCH
a = -60.65;
b = 0.0006635;
c = 38.16;
d = -0.06593;
nu_c = fs/2;

psd1side= biExpPsd(a,b,c,d,nu_c,freqAx);
pitchModelOut = getRotation(pitch, psd1side, Ns, df);

if isPlot
    subplot(3,1,2)
    plot((1:Ns)*4e-3, wrapTo180(pitchModelOut), '*')
end