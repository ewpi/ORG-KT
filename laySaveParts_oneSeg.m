% Code to generate binary files of layers 
% Segmentation for all BM-scans
% Code needs interaction with user

close all
clear all

load('Structure.mat');
rangeAscans = Structure.rangeAscans;

pix1 = 12;
howManyPixels = 150;

figures = 1;
savelayers = 1;

path01 = [cd '\Example\'];  % mouse, eye1s, Data
FilesList = dir([path01 '*Amp.bin']);
Nfiles = length(FilesList);

if contains(path01, 'mouse')
    howManyPixels = 200;
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
bscanpath = [bscanpath '.jpg'];

disp('Read data, amp file')

% *Read amplitude data* 
pathA = [path0 folder path3 path2];
ampc = read3DleeS(char(pathA));
ampc = single(ampc);
SA = size(ampc);

    aver_bscan = mean(ampc(:,:,:),3);

if figures == 1
    figure(10)
    clear image1
    image1 = aver_bscan;
    figure(); colormap(gray(255)); imagesc(20*log10(image1)); set(gca,'Ydir','normal');
    axis equal; colorbar; clim([40 90]); axis off   
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
BscanPreSegment(:,:) = squeeze(mean(ampc(:,:,S(3)/2-2:S(3)/2+2),3));
BscanPreSegment = 20*log10(BscanPreSegment);
BscanPreSegment = rescale(BscanPreSegment,0,255);
BscanPreSegment = uint8(BscanPreSegment); 
BscanPreSegment = flip(BscanPreSegment,1);
imshow(BscanPreSegment,[]);
imwrite(BscanPreSegment,bscanpath,'jpg');
%


% *Memory allocation*
array_phase = nan(howManyPixels,S(2),S(3));
array_amp = nan(howManyPixels,S(2),S(3));

array_phase = single(array_phase);
array_amp = single(array_amp);  

    % imPreSeg = imread(pathCSV);
    % %imshow(imPreSeg);
    % imPreSeg = rgb2gray(imPreSeg);
    % 
    % clear Kv
    % Kv = flip(sum(imPreSeg, 3), 1);
    % [~, presegv] = max(flip(Kv ~= 255, 1), [], 1);
    % presegv = size(Kv, 1) - presegv + 1;
    % presegv(~any(Kv ~= 255, 1)) = NaN;
    % preseg = median(presegv, 'omitnan');

    figure(11); imagesc(20*log10(aver_bscan)); colormap(gray(255));
        set(gca,'Ydir','normal'); clim([40 90]);
        title('Mark line ISOS, preliminary segmentation')
    
        % manual presegmentation
        h = drawpolyline('Color','y');
        pos = h.Position;
    
        presegv = interp1(pos(:,1),pos(:,2),1:S(2),'linear',NaN);
        presegv = round(presegv);
        preseg = median(presegv,'omitnan');

costv = nan(1,S(2));    
isosv = costv;


counter_frame = 0;
counter_bscan = 0;
for bscan = 1
    counter_bscan = counter_bscan + 1;

    counter_ascan = 0;
    for ascan = rangeAscans % ascans close to image borders should be cutted
    counter_ascan = counter_ascan + 1;    

    line1 = aver_bscan(:,ascan);

    line = movmean(line1, 2);
    halfWindow = round((pix1+1)/2.2);
    center = round(presegv(ascan));
    
    idx = max(1,center-halfWindow):min(length(line),center+halfWindow);
    
    [~, loc] = max(line(idx));
    x1w2_c = idx(loc) - 1;
    
    costv(ascan) = x1w2_c + pix1;
    isosv(ascan) = x1w2_c;

    end % ascans

    isosv_fit = polyval(polyfit(rangeAscans,isosv(rangeAscans),1),rangeAscans);
    %isosv_fit = reshape(feval(fit(rangeAscans(:), reshape(isosv(rangeAscans),[],1), 'exp1'), rangeAscans(:)), size(rangeAscans));
    %isosv_fit = reshape(feval(fit(rangeAscans(:), ...
    %reshape(isosv(rangeAscans),[],1), 'smoothingspline', 'SmoothingParam', 0.001), rangeAscans(:)), size(rangeAscans));
    costv_fit = isosv_fit+pix1;

     if figures == 1
        fig10 = figure('Visible','on');
        set(gcf,"Position",[10 50 1600 500])
        subplot(131); 
            plot((line), 'b-'); hold on; xline(center, 'y--', 'LineWidth', 2);
            xline(idx(1), 'k--'); xline(idx(end), 'k--'); locsM = idx(loc);
            plot(locsM, line(locsM), 'r+', 'LineWidth', 3, 'MarkerSize', 12); hold off;
            clear image2
            image2(:,:) = aver_bscan;
        subplot(132); 
            colormap(gray(255)); imagesc(20*log10(image1)); set(gca,'Ydir','normal');
            colorbar; clim([40 90]);
            hold on; plot(1:length(presegv),presegv,'y-','LineWidth',3,'MarkerSize',12);
            % axis equal
            plot(rangeAscans,isosv(rangeAscans),'r+','LineWidth',3,'MarkerSize',12);
            plot(rangeAscans,costv(rangeAscans),'b+','LineWidth',3,'MarkerSize',12);
            %plot(rangeAscans,isosv_fit,'b-','LineWidth',1.5);

        plot(rangeAscans,isosv(rangeAscans)-howManyPixels,'k+','LineWidth',3,'MarkerSize',12);
        plot(ascan,isosv(ascan),'g+','LineWidth',3,'MarkerSize',12); hold off

        subplot(133); 
            colormap(gray(255)); imagesc(20*log10(image1)); set(gca,'Ydir','normal');
            colorbar; clim([40 90]);
            % axis equal
            hold on; plot(rangeAscans,costv_fit,'b-','LineWidth',1.5);
            plot(rangeAscans,isosv_fit,'r-','LineWidth',1.5);

        xlim([ascan-50 ascan+50]); ylim([isosv_fit(1)-howManyPixels/2 isosv_fit(1)+howManyPixels/2])

        saveas(fig10,[pathSave 'SegLay_cost' num2str(round(howManyPixels/2)) '_isosDelta' num2str(pix1)])
     end

    disp('Segmentation done');    
end 

    counter_ascan = 0;
    array_phase = nan(howManyPixels,S(2),S(3),'single');
    array_amp = nan(howManyPixels,S(2),S(3),'single');

    for ascan = rangeAscans
        counter_ascan = find(rangeAscans == ascan);
    
        array_phase(:,ascan,:) = phasec(round(isosv_fit(counter_ascan))-howManyPixels/2:round(isosv_fit(counter_ascan))+howManyPixels/2-1,ascan,:);
        array_amp(:,ascan,:) = ampc(round(isosv_fit(counter_ascan))-howManyPixels/2:round(isosv_fit(counter_ascan))+howManyPixels/2-1,ascan,:);
        % cost layer howManyPixels/2, delta isos pix1
    end


if figures == 1
    figure(13)
        clear image1
        image1(:,:) = array_amp(:,100,:);
        figure(); colormap(gray(255)); imagesc(20*log10(image1)); 
        set(gca,'Ydir','normal');
        colorbar; clim([40 90]);  
end


disp('LaySave done');        
clear ampc phasec


if savelayers == 1
    path_array_phase = [pathSave path1 '_layers' path2];
    path_array_amp = [pathSave path3 '_layers' path2];

    save2DleeD(path_array_phase,double(array_phase));
    save2DleeD(path_array_amp,double(array_amp));
end


end