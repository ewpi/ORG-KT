close all
clear all

% Program to generate cross-spectra for A-scans and Pixels in time (in B-scans)
% Input data:
%           mainpath - main path to data and results
%           namedata - main filename of data 
%           pathData - file location
%           pathSaveB - path to saved cross-spectra
%           path_array_phase - path and filename to cut phase data with layers,
%                              saved from the first segmented layer with a constant distance from it
%           path_array_amp - path and filename to cut amplitude data with layers
%           lay_phase - phase of complex data saved from the first segmented layer with a constant distance from it,
%                       double precision, little-endian with size head in int32.
%                       Example data is lay_phase_example.bin
%           lay_amp - amplitude of complex OCT data saved from the first segmented layer with a constant distance from it,
%                       double precision, little-endian with size head in int32
%                       Example data is lay_amp_example.bin
%           read3DleeD - function to read lay_phase and lay_amp
%           save2DleeD - function to save cross-spectra


load('Structure.mat');
rangeAscans = Structure.rangeAscans;
lay2 = Structure.lay2;
selectLay2 = Structure.selectLay2;

range = [40 90];
path01 = [cd '\Example\'];  % mouse, eye1s

FilesList = dir([path01 '*Amp.bin']);
Nfiles = length(FilesList);

% OCT
for ind_file = 1:Nfiles
    
   filename = FilesList(ind_file).name;
   [filepath,base,ext] = fileparts(filename);

%folder = [base '_ZCor\'];
folder = [''];
pathSave = [path01 folder];

path11 = base;
path1 = strrep(base,'_Amp','_Phase');
namedata = strrep(base,'_Amp','');
path2 = '.bin';
pathData = path01;

pathSaveB = [path01 namedata '\lay\'];
if ~exist(pathSaveB, 'dir')
   mkdir(pathSaveB)
end

path_array_phase = [pathSave path1 '_layers' path2];
path_array_amp = [pathSave base '_layers' path2];


% The read range of data with interesting layers, 
% saved from the first segmented layer with a constant distance from it
%

lay_phase = read3DleeD(path_array_phase);
lay_amp = read3DleeD(path_array_amp); 
S3 = size(lay_phase,3);
S = size(lay_amp);

lineToSeg = squeeze(mean(lay_amp,[2,3],'omitnan'));
meanBscan = squeeze(mean(lay_amp,2,'omitnan'));
%lineToSeg(1:10)=0;
%lineToSeg(end-10:end)=0;
lineToSeg = lineToSeg/max(lineToSeg,[],"all");
[peaks,locs]=findpeaks(lineToSeg,'MinPeakHeight',0.38, 'MinPeakDistance',3);

flip_lineToSeg = flip(lineToSeg);
[peaksf,locsf]=findpeaks(flip_lineToSeg,'MinPeakHeight',0.38, 'MinPeakDistance',3);
locsf = S(1)-locsf;

lay1 = locs(1);
%lay2 = locsf(1);
% 
% lay2 = 100; %%%%% if segmentation does not work, put constant value, distans from ELM

Structure.rangeAscans=rangeAscans;
Structure.lay1 = lay1;
% Structure.lay2 = lay2;
Structure.S=S;

fig5 = figure();
subplot(211); colormap(gray(255)); imagesc(meanBscan);
subplot(212); plot(lineToSeg)
    hold on; plot(lay1,lineToSeg(lay1),'r+','LineWidth',3,'MarkerSize',12)
    ylim([0 1.2])

if selectLay2 == 1
    [xi,yi] = getpts(fig5);
    lay2 = round(xi);
    Structure.lay2 = lay2;
end
         plot(lay2,lineToSeg(lay2),'b+','LineWidth',3,'MarkerSize',12);
         subplot(211); colormap(gray(255)); imagesc(meanBscan);
             hold on; plot([1,S(3)],[lay2 lay2],'b-','LineWidth',1.5);

fig5 = figure();
subplot(211); colormap(gray(255)); imagesc(meanBscan);
subplot(212); plot(lineToSeg); ylim([0 1.2])

selectLay1 = 1;
if selectLay1 == 1
    [x2i,y2i] = getpts(fig5);
    lay1 = round(x2i);
    Structure.lay1 = lay1;
end   
         hold on; plot(lay1,lineToSeg(lay1),'r+','LineWidth',3,'MarkerSize',12)
         plot(lay2,lineToSeg(lay2),'b+','LineWidth',3,'MarkerSize',12);
         subplot(211); colormap(gray(255)); imagesc(meanBscan);
             hold on; plot([1,S(3)],[lay2 lay2],'b-','LineWidth',1.5);
             plot([1,S(3)],[lay1 lay1],'r-','LineWidth',1.5); hold off


saveas(gcf,[pathSave 'LayPosTest.jpg'],'jpg');
save([path01 namedata '\Structure.mat'],'Structure');

m = nan(S3,S3);
for r = rangeAscans % Range of interesting A-scans

    for p = [lay1,lay2] % Pixels indicate layers   
    % for p = [16,33] % Pixels indicate layers   
        m = nan(S3,S3);
        
        for bscan = 1:S3-1
             
            % Calculate columns of cross-spectra, 
            % the B-scan value is multiply by signal from this B-scan to end

            L = length(bscan:S3);
            w1 = lay_amp(p,r,bscan).*exp(1i*lay_phase(p,r,bscan));
            w2 = squeeze(lay_amp(p,r,bscan:S3).*exp(1i*lay_phase(p,r,bscan:S3)));
            
            m(1:L,bscan) = w1.*conj(w2(:));
            
            clear w1 w2 L
        end
        
        % ---- save cross-spectum for A-scan r and Pixel p----
            clear m_phase m_amp path_array_phase path_array_amp r1 mmax
            
        m_phase = angle(m);
        m_amp = abs(m);
        
        path_array_phase_cross = [pathSaveB namedata '_Phase' '_p' num2str(p) '_A' num2str(r)];
        path_array_amp_cross = [pathSaveB namedata '_Amp' '_p' num2str(p) '_A' num2str(r)];
        
        save2DleeD(path_array_phase_cross,m_phase);
        save2DleeD(path_array_amp_cross,m_amp); 
            
        close all
    p
    end
    r
end

end
