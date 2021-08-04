function angleModelOut = getRotation(angle, psd1side, Ns, df, varargin)
%GETROTATION returns the syntetic angle process.

%   GETROTATION(angle, psd1side, Ns, df, varargin) returns the time domain 
%   angle process as Section 3/Figure 5 in Steve Blandino, 
%   Tanguy Ropitault, Raied Caromi, Jacob Chakareski, Mahmudur Khan, 
%   and Nada Golmie. 2021. Head Rotation Model for VR System Level 
%   Simulations. 

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

if isempty(varargin)
    varNoise = 0;
else
    varNoise = varargin{1};
end
    
angle = angle(:);
awgn  = 1 + sqrt(varNoise) * (randn(Ns,1) );
fftAngle = fft(angle, 2*Ns);
fftAngle1side = fftAngle(1:end/2);
fftAngleFilter = (fftAngle1side.*sqrt(df*10.^(psd1side/10))/sqrt(df*10.^(psd1side(1)/10))).*awgn;%%/sqrt(df*10.^(psd1side(1)/10)))
angleModelOut = ifft([fftAngleFilter(1:end); conj(fftAngleFilter(end:-1:2))]);
angleModelOut = real(angleModelOut(1:Ns));

end