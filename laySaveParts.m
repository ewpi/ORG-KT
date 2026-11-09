% Code to generate binary files of layers 
% Segmentation for all BM-scans
% Code needs interaction with user

close all
clear all

load('Structure.mat');
rangeAscans = Structure.rangeAscans;

pix1 = 12;
howManyPixels = 150;

figures = 0;
savelayers = 1;

path01 = [cd '\Example\'];  % mouse, eye1s, Data
FilesList = dir([path01 '*Amp.bin']);
Nfiles = length(FilesList);

if contains(path01, 'mouse')
    howManyPixels = 120;
end

% OCT
for ind_file = 1:Nfiles
    
   filename = FilesList(ind_file).name;
    [filepath,base,ext] = fileparts(filename);

%folder = [base '_ZCor\'];
folder = [''];
pathSave = [path01 folder];

path0 = path01;
range = [40 90];

path11 = base;
path3 = base;
path1 = strrep(base,'_Amp','_Phase');
namedata = strrep(base,'_Amp','');
path2 = '.bin';

bscanpath = [path0 namedata '_PresegBscan'];
pathCSV = [bscanpath 'L.jpg'];
%pathCSV = [bscanpath 'L2.jpg'];
bscanpath = [bscanpath '.jpg'];

disp('Read data, amp file')

% *Read amplitude data* 
pathA = [path0 folder path3 path2];
ampc = read3DleeS(char(pathA));
ampc = single(ampc);
SA = size(ampc);

if figures == 1
    figure()
    clear image1
    image1(:,:) = ampc(:,:,1);
    figure(); colormap(gray(255)); imagesc(20*log10(image1)); set(gca,'Ydir','normal');
    axis equal; colorbar; caxis([40 90]); axis off
end
    
disp('Read data, phase file')

% *Read phase data* 
pathP = [path0 folder path1 path2];
phasec = read3DleeS(char(pathP));
phasec = single(phasec);
S = size(phasec);

% ampcTest(:,:,:) = ampc(:,:,1:2:S(3));
% phasecTest(:,:,:) = phasec(:,:,1:2:S(3));
% save2DleeS(pathA,ampcTest);
% save2DleeS(pathP,phasecTest);


% *Save B-scan to presegmentation (open jpg in Gimp and save layer)*
%
clear BscanPreSegment
BscanPreSegment(:,:) = squeeze(mean(ampc(:,:,S(3)/2-2:S(3)/2+2),3));
BscanPreSegment = 20*log10(BscanPreSegment);
BscanPreSegment = rescale(BscanPreSegment,0,255);
BscanPreSegment = uint8(BscanPreSegment); 
%BscanPreSegment = flip(BscanPreSegment,1);
imshow(BscanPreSegment,[]);
imwrite(BscanPreSegment,bscanpath,'jpg');
%


% *Memory allocation*
array_phase = nan(howManyPixels,S(2),S(3));
array_amp = nan(howManyPixels,S(2),S(3));

array_phase = single(array_phase);
array_amp = single(array_amp);  

    imPreSeg = imread(pathCSV);
    imshow(imPreSeg);
    imPreSeg = flipud(imcomplement(rgb2gray(imPreSeg)));

    figure()
    imagesc(imPreSeg); set(gca,'Ydir','normal');


K1_v = imPreSeg;
[~, preseg_v] = max(flipud(K1_v), [], 1);
preseg_v = double(size(K1_v, 1) - preseg_v + 1);
preseg_v(~any(K1_v, 1)) = NaN;

counter = 0;
counter_frame = 0;
counter_bscan = 0;
for bscan = 1:1:S(3)
    counter_bscan = counter_bscan + 1;

    for ascan = rangeAscans % ascans close to image borders should be cutted
    counter = counter + 1;    

    preseg = preseg_v(ascan);
    
    if bscan > S(3)-3
        line1 = squeeze(mean(ampc(:,ascan-3:ascan+3,bscan-2:end),[2,3]));  
    elseif bscan-3 < 0
        line1 = squeeze(mean(ampc(:,ascan-3:ascan+3,bscan:bscan+2),[2,3])); 
    else
        line1 = squeeze(mean(ampc(:,ascan-3:ascan+3,bscan-2:bscan+2),[2,3]));   
    end


    windowSize = 3; 
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;    
    line = filter(b,a,line1);    
    line(1:20) = 0;
    line(end-20) = 0;
    line(preseg+round(pix1/2):end) = 0;
    line(1:preseg+-round(pix1/2)) = 0;


    if figures == 1
        fig10 = figure('Visible','on');
        set(gcf,"Position",[10 50 1600 500])
    end

    line = line/max(line,[],'all');
    line2 = flip(line);
    [peaks,locs]=findpeaks(line2,'MinPeakHeight',0.35, 'MinPeakDistance',3);
    locs = S(1)-locs;
    [Mmax, locsM] = max(line,[],'all');
    x1w2_c = locsM;
        
    if figures == 1
        
        subplot(131); plot(line,'y-');    
        hold on; subplot(131); plot(line); hold off
        hold on; plot(line2,'g-');
        hold on; plot(locsM,line(locsM),'g+','LineWidth',3,'MarkerSize',12); hold off
        hold on; plot(x1w2_c,line(x1w2_c),'r+','LineWidth',3,'MarkerSize',12); hold off

        clear image1
        image1(:,:) = ampc(:,:,bscan);
        subplot(132); colormap(gray(255)); imagesc(20*log10(image1)); set(gca,'Ydir','normal');
        colorbar; caxis([40 90]);
        % axis equal
        hold on; plot(ascan,x1w2_c,'r+','LineWidth',3,'MarkerSize',12); hold off
        hold on; plot(ascan,x1w2_c+pix1,'b+','LineWidth',3,'MarkerSize',12); hold off
        hold on; plot(ascan,x1w2_c+pix1-howManyPixels,'k+','LineWidth',3,'MarkerSize',12); hold off
        hold on; plot(ascan,locsM,'g+','LineWidth',3,'MarkerSize',12); hold off

        hold on; plot(1:length(preseg_v),preseg_v,'y-','LineWidth',1.5); hold off

        subplot(133); colormap(gray(255)); imagesc(20*log10(image1)); set(gca,'Ydir','normal');
        colorbar; caxis([40 90]);
        % axis equal
        hold on; plot(ascan,x1w2_c,'r+','LineWidth',3,'MarkerSize',12); hold off
        hold on; plot(ascan,x1w2_c+pix1,'b+','LineWidth',3,'MarkerSize',12); hold off
        hold on; plot(ascan,x1w2_c+pix1-howManyPixels,'k+','LineWidth',3,'MarkerSize',12); hold off

        xlim([ascan-50 ascan+50]); ylim([x1w2_c-howManyPixels/2 x1w2_c+howManyPixels/2])

        if mod(counter,10) == 0
            counter_frame = counter_frame + 1;
            title(['Bscan ' num2str(bscan) ', Ascan ' num2str(ascan)]);
            Frame1(counter_frame) = getframe(fig10);
        end
        close all
    end     
        
        array_phase(1:howManyPixels,ascan,counter_bscan) = phasec((x1w2_c-howManyPixels/2):(x1w2_c+howManyPixels/2)-1,ascan,bscan);
        array_amp(1:howManyPixels,ascan,counter_bscan) = ampc((x1w2_c-howManyPixels/2):(x1w2_c+howManyPixels/2)-1,ascan,bscan);
       
    end
    bscan
end 
clear ampc phasec

end

if figures == 1
    figure(13)
        clear image1
        image1(:,:) = array_amp(:,100,:);
        figure(); colormap(gray(255)); imagesc(20*log10(image1)); 
        set(gca,'Ydir','normal');
        colorbar; clim([40 90]);  
end

if savelayers == 1
    path_array_phase = [pathSave path1 '_layers' path2];
    path_array_amp = [pathSave path3 '_layers' path2];

    save2DleeD(path_array_phase,double(array_phase));
    save2DleeD(path_array_amp,double(array_amp));
end

if figures == 1
    v = VideoWriter([pathSave namedata '_segTest' '.avi']);
    v.FrameRate = 4;
    open(v);
    writeVideo(v,Frame1)
    close(v);
end


