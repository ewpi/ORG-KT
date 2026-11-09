close all
clear all

names = {'183915_10s'};

S3 = 512;
NN = 20;
weights_binary_Yes_No = 1;
phaseRange1 = 2*pi;
phaseRange2 =10;  
phaseRange3 = 10;
    th = 0.16;  % threshold from weights
weightsKw = 1;
drawCS = 1;
mean_gaps = 1;
    pixelsToAv = 6;
freq = 50;
tmax = S3/freq;
dt = linspace(1,tmax,S3);
delta_t = dt(2)-dt(1);
m = nan(S3,S3);
pix1 = 43; 
pix2 = 99;
column = 1;
xx = linspace(0,tmax,S3);

for j=1
clear rangeAscans   
    if j == 1
    rangeAscans = 300:360;
    end
RA = length(rangeAscans);

for num_paths = 3

scaleToPosition0 = (1/860e-9)*(4*1.333333*pi); % [1/m]
scaleToPosition = (1/(860))*(4*1.333333*pi); % mm
scaleToPosition1 = (1/(860e-9)/1e6)*(4*1.333333*pi); % mm

numfile = 0;
for i=names
    
numfile = numfile + 1;    
name0 = char(i);

PathSigma = cd;
MainPath = [PathSigma '\Example\' name0 '\lay12\'];    
pathFExt = '.bin';
pathData = [MainPath]; 

pathSaveB = [PathSigma '\Example\' name0 '\org\'];
if ~exist(pathSaveB, 'dir')
   mkdir(pathSaveB)
end

name = [name0];

path_phase = [MainPath name 'Phase_lay1_' num2str(pix1) '_lay2_' num2str(pix2)...
             '_RA' num2str(RA) '_fromAsc_' num2str(rangeAscans(1)) '_lay12dif.bin'];
path_w = [MainPath name '_weights_lay1_' num2str(pix1) '_lay2_' num2str(pix2)...
         '_RA' num2str(RA) '_fromAsc_' num2str(rangeAscans(1)) '_layDif.bin'];
m_phase = read2DleeD(path_phase);
aw = read2DleeD(path_w);

m_phase = m_phase(1:S3,1:S3);
aw = aw(1:S3,1:S3);

A = ones(size(m_phase));
A = tril(A);
A = flip(A,1);
kk = find(A==0);
B = m_phase.*A;
B(kk)=nan;
m_phase = B;
aw = aw.*A;
aw(kk)=nan;

% ************ FIGURES ************
C = m_phase(:,1);
r1 = find(isnan(m_phase)); 
r2 = find(isnan(m_phase)==0);
m_phase(r1) = 0;

aw_raw = aw.*aw;
if weightsKw == 1
    aw = aw.*aw;
end
% binary weights
if weights_binary_Yes_No == 1
    ind2 = find(aw<th);
    ind3 = find(aw>th);
    aw(ind2) = 0;
    aw(ind3) = 1;
elseif weights_binary_Yes_No == 3   
    ind2 = find(aw<0.40 & aw>0.1);
    aw(ind2) = 0.15;
end

m_phase_aw = m_phase.*aw;


% STD first column CD
[c,edge] = histcounts(C,length(-pi:0.01:pi));
pro_z = c;
fi = edge(1:end-1);
    % linear std
    mas_o = sum(pro_z);
    mea_o = sum(pro_z.*fi)/mas_o;
    var_o = sum(pro_z.*(fi-mea_o).^2)/mas_o;
    sig_o = sqrt(var_o);
    sig_o_pos1 = sig_o/scaleToPosition;

    % circular std
    sum_c = sum(pro_z.*exp(1i*fi))/mas_o;
    sig_c1 = sqrt(-2*log(abs(sum_c)));
    sig_c1_pos = sig_c1/scaleToPosition;
    
    if numfile > 5
        correction2 = sig_c1/4;
    else
        correction2 = 0;
    end

% ----------- PLOT PHA 1 --------------------
    s3 = m_phase(:,column);
    f2 = figure('Visible','on'); 
    yyaxis left
    p1 = plot(xx,m_phase(:,column),'b-','LineWidth',2);
    KT = ['Knox-Thompson paths ' num2str(num_paths)];
    KT2 = ['KT paths ' num2str(num_paths)];
    hold on; plot(xx,zeros(size(xx)),'k--'); hold off
    hold on; plot(xx,pi*ones(size(xx)),'k--'); hold off
    hold on; plot(xx,-pi*ones(size(xx)),'k--'); hold off
    legend([p1],[num2str(column) ' (st) column of arg(<U_{12}>)'],...
          'Location','northoutside','Orientation','horizontal')
    legend boxoff; ylim([-phaseRange1 phaseRange2])
    set(f2,'Position',[200 150 1200 400]); xlim([xx(1) xx(end)]);
    ylabel('Phase \phi'), xlabel('Time t [s]');
    set(gca,'TickLabelInterpreter','latex','FontSize', 20, 'FontName', 'Times New Roman')
    yyaxis right
    ylim([-phaseRange1/scaleToPosition phaseRange2/scaleToPosition]); ylabel('Shift [nm]')
    text(8.3,400,['\sigma_{CD} = ' num2str(sig_c1_pos,'%.02f')],'FontSize',20,'FontName','Times New Roman','Color',[1,0.37,0])
    set(gca,'TickLabelInterpreter','latex','FontSize', 20,'FontName','Times New Roman')
       saveas(f2,char([pathSaveB 'w1_KT_CD_Plot_paths' num2str(num_paths) '_RA' num2str(RA) '_pix1_' num2str(pix1)...
       '_pix2_' num2str(pix2) '_fromAsc_' num2str(rangeAscans(1)) '.jpg']),'jpg')


% Knox-Thomson calculatons
del_paths = floor(S3/num_paths);
% if mod(S3,num_paths) == 0
    colap1 = m_phase(1:del_paths,1);
    wcolap1 = m_phase_aw(1:del_paths,1);
    line_aw =  aw(1:del_paths,1);
% end   
for il = 1:num_paths-1

        corrector = mean(colap1(end-NN:end))+correction2;
        correctorw = mean(wcolap1(end-NN:end))+correction2;

    colap_next = m_phase(1:del_paths,il*del_paths)+corrector;
    colap1 = [colap1; colap_next];
    wcolap_next = m_phase_aw(1:del_paths,il*del_paths)+correctorw;
    wcolap1 = [wcolap1; wcolap_next];
    w_next = aw(1:del_paths,il*del_paths);
    line_aw  = [line_aw; w_next];
end

% STD from Knox-Thomson paths
KT1 = colap1./scaleToPosition;
clear c edge
[c,edge] = histcounts(colap1,length(-pi:0.01:pi));
pro_z = c;
fi = edge(1:end-1);
   %linear z momentów
    mas_o = sum(pro_z);
    mea_o = sum(pro_z.*fi)/mas_o;
    var_o = sum(pro_z.*(fi-mea_o).^2)/mas_o;
    sig_o = sqrt(var_o);
    sig_o_pos = sig_o/scaleToPosition;

    %circular 
    sum_c = sum(pro_z.*exp(1i*fi))/mas_o;
    mea_c = angle(sum_c);
    var_c = 1-abs(sum_c);
    sig_c = sqrt(-2*log(abs(sum_c)));
    sig_c_pos = sig_c/scaleToPosition;


% ----------- PLOT PHA 2 --------------------
    f2 = figure('Visible','on'); 
    p1 = plot(xx(1:length(colap1)),colap1,'r-','LineWidth',2);
    KT = ['Knox-Thompson paths ' num2str(num_paths)];
    wKT = ['weight and KT paths ' num2str(num_paths)];
    hold on; plot(xx,zeros(size(xx)),'k--'); hold off
    hold on; plot(xx,pi*ones(size(xx)),'k--'); hold off
    hold on; plot(xx,-pi*ones(size(xx)),'k--'); hold off
    hold on; p4 = plot(xx(1:length(colap1)),colap1.*line_aw,'k-','LineWidth',1); hold off
    legend([p1, p4],[KT],['Weights\cdot(' KT2 ')'],...
          'Location','northoutside','Orientation','horizontal')
    legend boxoff; ylim([-phaseRange1 phaseRange2])
    set(f2,'Position',[200 150 1200 400]); xlim([xx(1) xx(end)]);
    ylabel('Phase \phi'), xlabel('Time t [s]');
    set(gca,'TickLabelInterpreter','latex','FontSize', 20,'FontName', 'Times New Roman')
       saveas(f2,char([pathSaveB 'w2_KT_CD_Plot_paths' num2str(num_paths) '_RA' num2str(RA) '_pix1_' num2str(pix1)...
       '_pix2_' num2str(pix2) '_fromAsc_' num2str(rangeAscans(1)) '.jpg']),'jpg')


if mean_gaps == 1
    gaps = find(line_aw == 0);
    vInGaps = colap1;
    N = length(gaps);
    for ii=1:N
        v11 = gaps(ii)-pixelsToAv;
        v22 = gaps(ii)+pixelsToAv;
        if v11 < 1
            v11 = gaps(ii);
        end
        if v22 > length(colap1)
            v22 = gaps(ii);
        end
        vInGaps(gaps(ii))=mean(colap1(v11:v22));
    end
    
    
% STD from Knox-Thomson paths and gaps
clear c edge sum_c mea_c var_c pro_z
[c,edge] = histcounts(vInGaps,length(-pi:0.01:pi));
pro_z = c;
fi = edge(1:end-1);
    % linear std
    mas_o = sum(pro_z);

    % circular std
    sum_c = sum(pro_z.*exp(1i*fi))/mas_o;
    mea_c = angle(sum_c);
    var_c = 1-abs(sum_c);
    sig_c3 = sqrt(-2*log(abs(sum_c)));
    sig_c3_pos = sig_c3/scaleToPosition;    
    
    % ----------- PLOT PHA 2 --------------------

    % ##
    f22 = figure('Visible','on'); 
    yyaxis left
    p1 = plot(xx(1:length(colap1)),colap1,'r-','LineWidth',2);
    hold on; plot(xx,zeros(size(xx)),'k--'); hold off
    hold on; plot(xx,pi*ones(size(xx)),'k--'); hold off
    hold on; plot(xx,-pi*ones(size(xx)),'k--'); hold off
    hold on; p4 = plot(xx(1:length(colap1)),vInGaps,'k-','LineWidth',1); hold off
    leg22 = legend([p1, p4],KT2,'Weighted ORG',...
          'Location','northoutside','Orientation','horizontal');
    legend boxoff; ylim([-phaseRange1 phaseRange3])
    set(leg22,'FontSize',20);
    set(f22,'Position',[200 150 600 400]); xlim([xx(1) 10]);
    ylabel('Phase \phi'), xlabel('Time t [s]');
%     xtickformat('%.01f')

    yyaxis right
    % hold on; p5 = plot(xx(1:length(vInGaps)),vInGaps/scaleToPosition,'k-','LineWidth',1); hold off
    ylim([-phaseRange1/scaleToPosition phaseRange3/scaleToPosition]); ylabel('Shift [nm]')
    ytickformat('%i')
    %text(4.3,0.24,['\sigma_{KT} = ' num2str(sig_c_pos,'%.02f')],'FontSize',20,'FontName','Times New Roman','Color',[1,0.37,0])

    set(gca,'TickLabelInterpreter','latex','FontSize', 20,'FontName','Times New Roman')
       saveas(f22,char([pathSaveB 'g22_KT_CD_Plot_paths' num2str(num_paths) '_RA' num2str(RA) '_pix1_' num2str(pix1)...
       '_pix2_' num2str(pix2) '_fromAsc_' num2str(rangeAscans(1)) '.jpg']),'jpg')
       saveas(f22,char([pathSaveB 'g22_KT_CD_Plot_paths' num2str(num_paths) '_RA' num2str(RA) '_pix1_' num2str(pix1)...
       '_pix2_' num2str(pix2) '_fromAsc_' num2str(rangeAscans(1)) '.fig']),'fig')

  
   % ##
    f22 = figure('Visible','on'); 
    yyaxis left
    
    KT = ['Knox-Thompson paths ' num2str(num_paths)];
    wKT = ['weight and KT paths ' num2str(num_paths)];
    plot(xx,zeros(size(xx)),'k--');
    hold on; plot(xx,pi*ones(size(xx)),'k--'); hold off
    hold on; plot(xx,-pi*ones(size(xx)),'k--'); hold off
    hold on; p4 = plot(xx(1:length(colap1)),vInGaps,'k-','LineWidth',1); hold off
    leg22 = legend([p4],'Weighted ORG',...
          'Location','northoutside','Orientation','horizontal');
    legend boxoff; ylim([-phaseRange1 phaseRange3])
    set(leg22,'FontSize',20);
    set(f22,'Position',[200 150 600 400]); xlim([xx(1) 10]);
    ylabel('Phase \phi'), xlabel('Time t [s]');
    % xtickformat('%.01f')

    yyaxis right
    % hold on; p5 = plot(xx(1:length(vInGaps)),vInGaps/scaleToPosition,'k-','LineWidth',1); hold off
    ylim([-phaseRange1/scaleToPosition phaseRange3/scaleToPosition]); ylabel('Shift [nm]')
 %    ytickformat('%.01f')
    %text(4.3,0.24,['\sigma_{KT} = ' num2str(sig_c_pos,'%.02f')],'FontSize',20,'FontName','Times New Roman','Color',[1,0.37,0])

    set(gca,'TickLabelInterpreter','latex','FontSize', 20,'FontName','Times New Roman')
       saveas(f22,char([pathSaveB 'g33_KT_CD_Plot_paths' num2str(num_paths) '_RA' num2str(RA) '_pix1_' num2str(pix1)...
       '_pix2_' num2str(pix2) '_fromAsc_' num2str(rangeAscans(1)) '.jpg']),'jpg')
       saveas(f22,char([pathSaveB 'g33_KT_CD_Plot_paths' num2str(num_paths) '_RA' num2str(RA) '_pix1_' num2str(pix1)...
       '_pix2_' num2str(pix2) '_fromAsc_' num2str(rangeAscans(1)) '.fig']),'fig')

end


if drawCS == 1
    % ----------- PLOT PHA 1 ------------------
    clear im; im = m_phase; 
    im(r1) = 0;
    f1 = figure('Visible','off'); colormap(mapaRWB); imagesc(dt,dt,im); colorbar; 
    caxis([-pi pi]);
    % caxis([-0.3 0.3]);
    set(f1,'Position',[200 150 1000 800]);
    set(gca,'Ydir','normal');
    set(gca,'TickLabelInterpreter','latex','FontSize', 32,'FontName','Times New Roman')
    xlabel('t'); ylabel('\Deltat');
    saveas(f1,char([pathSaveB 'P_paths' num2str(num_paths) '_RA' num2str(RA) '_pix1_' num2str(pix1)...
    '_pix2_' num2str(pix2) '_fromAsc_' num2str(rangeAscans(1)) '.jpg']),'jpg')
    %saveas(f1,char([pathSaveB 'P' '.fig']),'fig')

    % ----------- PLOT PHA UNW ------------------
    clear im; im = m_phase_aw; 
    im(r1) = 0;
    f10 = figure('Visible','off'); colormap(mapaRWB); imagesc(dt,dt,im); colorbar; 
    caxis([-pi pi]);
    % caxis([-0.3 0.3]);
    set(f10,'Position',[200 150 1000 800]);
    set(gca,'Ydir','normal');
    set(gca,'TickLabelInterpreter','latex','FontSize', 32,'FontName','Times New Roman')
    xlabel('t'); ylabel('\Deltat');
    %saveas(f10,char([pathSaveB 'WP' '.fig']),'fig')    
    saveas(f10,char([pathSaveB 'WP_paths' num2str(num_paths) '_RA' num2str(RA) '_pix1_' num2str(pix1)...
       '_pix2_' num2str(pix2) '_fromAsc_' num2str(rangeAscans(1)) '.jpg']),'jpg')  

    % ----------- PLOT AMP 1 ------------------
        
        clear im; im = aw; 
        im(r1) = 1;
        f3 = figure('visible','on'); colormap(gray(2)); imagesc(dt,dt,im); col = colorbar; 
        set(col,'XTick',[0,1]);
    %     col.Ruler.TickLabelFormat = '%.1f';
    %     caxis([0 1e14]);
        set(f3,'Position',[200 150 1000 800]);
        set(gca,'Ydir','normal');
        set(gca,'TickLabelInterpreter','latex','FontSize', 32,'FontName','Times New Roman')
        xlabel('t'); ylabel('\Deltat');
        saveas(f3,char([pathSaveB 'Wbinarne_paths' num2str(num_paths) '_RA' num2str(RA) '_pix1_' num2str(pix1)...
       '_pix2_' num2str(pix2) '_fromAsc_' num2str(rangeAscans(1)) '.jpg']),'jpg')
        %saveas(f3,char([pathSaveB 'Wbinarne' '.fig']),'fig')

        clear im; im = aw_raw; 
        im(r1) = 1;
        f3 = figure('visible','off'); colormap(hot); imagesc(dt,dt,im); col = colorbar; 
   %     col.Ruler.TickLabelFormat = '%.1f';
        set(f3,'Position',[200 150 1000 800]);
        set(gca,'Ydir','normal');
        set(gca,'TickLabelInterpreter','latex','FontSize', 32,'FontName','Times New Roman')
        xlabel('t'); ylabel('\Deltat');
        saveas(f3,char([pathSaveB 'Wniebinarne_paths' num2str(num_paths) '_RA' num2str(RA) '_pix1_' num2str(pix1)...
       '_pix2_' num2str(pix2) '_fromAsc_' num2str(rangeAscans(1)) '.jpg']),'jpg')
        %saveas(f3,char([pathSaveB 'Wniebinarne' '.fig']),'fig')
        
 % how many percentage we cut
 aw(r1)=nan;
 pozostaleAll = find(aw==1);
 cpozos = length(pozostaleAll);
 allP = length(r2);
 procentAll = cpozos/allP*100
 pozostaleOne = find(line_aw);
 oneP = length(pozostaleOne);
 procentOne = oneP/S3*100
 
 
    f44 = figure(44);
    fx = linspace(1/tmax,freq,length(colap1)); 
    set(f44,'Position',[100 100 600 400]);
    s1 = abs(fft(exp(1i*colap1)));
    loglog(fx,s1,'r-');
    ylabel('|FFT(exp(\Delta\phi))| [dB]'); xlim([fx(1) fx(end)])
    xlabel('Frequency [Hz]');
    set(gca,'TickLabelInterpreter','latex','FontSize', 20,'FontName','Times New Roman')
    saveas(f44,char([pathSaveB 'FFT_paths' num2str(num_paths) '_RA' num2str(RA) '_pix1_' num2str(pix1)...
       '_pix2_' num2str(pix2) '_fromAsc_' num2str(rangeAscans(1)) '.jpg']),'jpg')
   
   
    f55 = figure(44);
    fx = linspace(1/tmax,freq,length(colap1)); 
    set(f55,'Position',[100 100 600 400]);
    s2 = abs(fft(exp(1i*vInGaps)));
    
    [f_max,f_max_pos] = max(s2(100:170));
    [sf_max,sf_max_pos] = max(s2(171:200));
    sf_freq(2,numfile)=fx(sf_max_pos+170);
    
    f_freq(2,numfile)=fx(f_max_pos+99);
    s22 = fft(exp(1i*vInGaps));
    s22(f_max_pos+99)=0;
    is22 = angle(ifft(s22));
    
    s22(sf_max_pos+170-3:sf_max_pos+170+3)=0;
    is33 = angle(ifft(s22));
    
    loglog(fx,s2,'k-');
    %hold on; loglog(fx(f_max_pos+99),f_max,'b+','MarkerSize',15); hold off
    ylabel('|FFT(exp(\Delta\phi))| [dB]'); xlim([fx(1) fx(end)])
    xlabel('Frequency [Hz]');
    set(gca,'TickLabelInterpreter','latex','FontSize', 20,'FontName','Times New Roman')
    saveas(f55,char([pathSaveB 'FFT2_paths' num2str(num_paths) '_RA' num2str(RA) '_pix1_' num2str(pix1)...
       '_pix2_' num2str(pix2) '_fromAsc_' num2str(rangeAscans(1)) '.jpg']),'jpg')
   
    ts1(numfile,:) = s1;
    ts2(numfile,:) = s2;
    org(numfile,:) = vInGaps; 
    ts3(numfile,:) = s3;
    ts11(numfile,:) = is22;
    ts33(numfile,:) = is33;
end
end
   
    s = vInGaps(200:end)/scaleToPosition;
    vx = xx(1:length(s));    
    clear P
    P = polyfit(vx,s',1);
    figure()
    plot(vx,s)
    hold on; plot(vx,P(1)*vx+P(2)); hold off
    v = s';
    diff = v-(P(1).*vx + P(2));
    hold on; plot(vx,diff); hold off
    hold on; plot(vx,zeros(size(vx))); hold off
    ylabel('[\mum]')
    resultSTD = std(diff)
    title(['STD = ' num2str(resultSTD,'%.1f') ' [nm]'])

    F_diff = (fft(diff));
    F_diff_a = abs(F_diff);

    F_diff_2 = (fft(s));
    mask_F_diff_2 = ones(size(F_diff_2));
    S_F_diff_2 = length(F_diff_2);
    ll = 1; pp = 10;
    mask_F_diff_2(ll:pp)=0;
    mask_F_diff_2(S_F_diff_2-pp:S_F_diff_2-ll)=0;

%     ll = 125; pp = 135;
%     F_diff_2(ll:pp)=0;
%     F_diff_2(S_F_diff_2-pp:S_F_diff_2-ll)=0;


    mask_F_diff_2 = mask_F_diff_2;
    F_diff_2 = F_diff_2.*mask_F_diff_2;


    figure();
    subplot(221); plot(s)
    subplot(222); plot(F_diff_a)
    subplot(223); plot(abs(F_diff_2))
    subplot(224); plot(abs(ifft(F_diff_2)))

    
end

clear ts1 s1 ts2 s2 org vInGaps ts3 s3 ts11 is22

end

