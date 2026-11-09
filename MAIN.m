% **** version MATLAB R2023b ***
% INPUT
%
% 1. Put your OCT processed data into Example folder.
% Binary files, *_Phase.bin and *_Amp.bin (3xint32 with size, and data in single/float numbers)
% Run the MAIN code.
% 2. Program needs interactions with user. In folder Example, .../Example/filename_PresegBscan.jpg
% you can find *jpg file. Open it in Gimp and create presegmentation of
% COST. Save your result at the same size in jpg format, .../Example/filename_PresegBscanL.jpg
% 3. run the MAIN code again


clear all
close all

% INPUT
rangeAscans = 100:120;    % A-scan range ( <50, RAM limit)
freq = 500;               % BM-scan frequency
lambda = 810e-9;          % lambda to scale phase -> nm
KT_paths = 2;             % Knox-Thompson paths
selectLay2 = 1;           % 1 - put second layer by hand on profile and press enter, 0 - use lay2 as a constant
    lay2 = 26;            % second layer, from presegmented retinal layer (ELM)
                          % blue marker on the LayPosTest.jpg
%
% end INPUT


Structure.rangeAscans=rangeAscans;
Structure.freq = freq;
Structure.lambda = lambda;
Structure.KT_paths = KT_paths;
Structure.lay2 = lay2;
Structure.selectLay2 = selectLay2;

save('Structure.mat','Structure');

cd
%run('laySaveParts_oneSeg.m')
run('laySaveParts.m')
run('cross_spectra_generator.m')
run('averaged_cross_spectra.m')


% OUTPUT
% cross-spectra phase average of 50th A-scans, COST and BrM

