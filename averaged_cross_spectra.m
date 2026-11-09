close all
clear all

% Program to calculate average of cross-spectra difference of layers     

figures = 0;
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
pathData = [path01 namedata '\lay\'];

load([path01 namedata '\Structure.mat']);

Bscans = Structure.S(3);
rangeAscans = Structure.rangeAscans;
lay1 = Structure.lay1;
lay2 = Structure.lay2;
freq = Structure.freq;
KT_paths = Structure.KT_paths;

total_time = Bscans/freq;
column = 1;
fx = freq/Bscans:freq/Bscans:freq;
xx = linspace(0,total_time,Bscans);
RA = length(rangeAscans);

pathSaveB = [path01 namedata '\lay12\'];
if ~exist(pathSaveB, 'dir')
   mkdir(pathSaveB)
end


licz = 0;
for ascan = rangeAscans
    licz = licz + 1;
        
    path_array_phase1 = [pathData namedata '_Phase_p' num2str(lay1) '_A' num2str(ascan)];
    path_array_amp1 = [pathData namedata '_Amp_p' num2str(lay1) '_A' num2str(ascan)];
    path_array_phase2 = [pathData namedata '_Phase_p' num2str(lay2) '_A' num2str(ascan)];
    path_array_amp2 = [pathData namedata '_Amp_p' num2str(lay2) '_A' num2str(ascan)];
    
    
    P1(:,:) = read2DleeD(path_array_phase1);
    A1(:,:) = read2DleeD(path_array_amp1);
    
    P2(:,:) = read2DleeD(path_array_phase2);
    A2(:,:) = read2DleeD(path_array_amp2);
    
    signalc(:,:,licz) = conj(A2.*exp(1i.*P2)) .* (A1.*exp(1i.*P1));
end


% =========== WEIGHTS ========================================

r1 = find(isnan(P1(:,:))); 
clear A1 P2 A2 m

w(:,:) = abs(sum(signalc,3))./sum(abs(signalc),3);
rnb = find(w<0);
w(r1)=0;
aw = w.^2;
clear w

path_w = [pathSaveB namedata '_weights_lay1_' num2str(lay1) '_lay2_' num2str(lay2) '_RA' num2str(RA) ...
'_fromAsc_' num2str(rangeAscans(1)) '_layDif.bin'];
save2DleeD(path_w,aw);
    
% =========== AVERAGE OF CROSS-SPECTRA and LAY DIFF ===========

signal = mean(signalc,3);
clear signalc
m_phase = angle(signal);
m_amp = abs(signal);

path_phase = [pathSaveB namedata 'Phase_lay1_' num2str(lay1) '_lay2_' num2str(lay2) '_RA' num2str(RA)...
'_fromAsc_' num2str(rangeAscans(1)) '_lay12dif.bin'];
path_amp = [pathSaveB namedata 'Amp_lay1_' num2str(lay1) '_lay2_' num2str(lay2) '_RA' num2str(RA)...
'_fromAsc_' num2str(rangeAscans(1)) '_lay12dif.bin'];

save2DleeD(path_phase,m_phase);
save2DleeD(path_amp,m_amp);
    
    % ************ FIGURES *****************************************
    
    m_phase(r1) = 0;
    aw(r1) = max(max(aw));
    mmax = max(max(m_amp));
    m_amp(r1) = mmax+10; 
    
    % ----------- PLOT PHASE ------------------
    clear path_array_phase path_array_amp
    
    f1 = figure('Visible','on'); colormap(mapaRWB); imagesc(xx,xx,m_phase); colorbar; 
    caxis([-pi pi]);
    set(f1,'Position',[200 10 1000 800]);
    set(gca,'Ydir','normal');
    set(gca,'TickLabelInterpreter','latex','FontSize', 28,'FontName','Times New Roman')
    xlabel('t [s]'); ylabel('\Deltat');
    saveas(f1,char([pathSaveB namedata 'Phase_lay1_' num2str(lay1) '_lay2_' num2str(lay2) '_RA' num2str(RA)...
    '_fromAsc_' num2str(rangeAscans(1)) '.jpg']),'jpg')
    
if figures == 1           
    % ----------- PLOT Phase*Weights ------------------
    m_phase_aw = m_phase.*aw;
    m_phase_aw(r1) = 0;    
    f10 = figure('Visible','on'); colormap(mapaRWB); imagesc(xx,xx,m_phase_aw); colorbar; 
    caxis([-pi pi]);
    set(f10,'Position',[200 10 1000 800]);
    set(gca,'Ydir','normal');
    set(gca,'TickLabelInterpreter','latex','FontSize', 28,'FontName','Times New Roman')
    xlabel('t [s]'); ylabel('\Deltat');
        saveas(f10,char([pathSaveB namedata 'wPhase_lay1_' num2str(lay1) '_lay2_' num2str(lay2) '_RA' num2str(RA)...
                   '_fromAsc_' num2str(rangeAscans(1)) '.jpg']),'jpg')    
           
    
    % ----------- WYKRES AMPLITUDE ------------------
    f3 = figure('visible','on'); colormap(hot); imagesc(xx,xx,m_amp); col = colorbar; 
    col.Ruler.TickLabelFormat = '%.1f';
    %caxis([0 1e14]);
    set(f3,'Position',[200 10 1000 800]);
    set(gca,'Ydir','normal');
    set(gca,'TickLabelInterpreter','latex','FontSize', 28,'FontName','Times New Roman')
    xlabel('t [s]'); ylabel('\Deltat');
       saveas(f3,char([pathSaveB namedata 'Amp_lay1_' num2str(lay1) '_lay2_' num2str(lay2) '_RA' num2str(RA)...
                   '_fromAsc_' num2str(rangeAscans(1)) '.jpg']),'jpg')
       
end    
    % % ----------- PLOT PHASE --------------------
    % f2 = figure('Visible','on','Color','w');
    % 
    % p1 = plot(xx,m_phase(:,column),'b-','LineWidth',1.8);
    % 
    % del_podz = floor(Bscans/KT_paths);
    % colap1 = m_phase(1:del_podz,1);
    % 
    % for il = 1:KT_paths-1
    %     colap_next = m_phase(1:del_podz,il*del_podz);
    %     colap1 = [colap1; colap_next + colap1(end)];
    % end
    % 
    % hold on;
    % p2 = plot(xx(1:length(colap1)),colap1,'r-','LineWidth',1.8);
    % hold off;
    % 
    % legend([p1,p2],'Direct calculation (reference t_0 = 0)', ...
    %     sprintf('Knox-Thompson (%d paths)',KT_paths), ...
    %     'Location','northoutside','Orientation','horizontal', ...
    %     'Box','off','FontName','Times New Roman','FontSize',15);
    % 
    % set(f2,'Position',[200 150 1100 420]);
    % set(gca,'FontName','Times New Roman','FontSize',17, ...
    %     'LineWidth',1.2,'Box','off','TickDir','out');
    % 
    % grid on;
    % xlim([xx(1) xx(end)]);
    % 
    % xlabel('Time, t (s)','FontName','Times New Roman','FontSize',19);
    % ylabel('Phase difference, \Delta\phi (rad)', ...
    %     'FontName','Times New Roman','FontSize',19);
    % 
    % exportgraphics(f2,[pathSaveB 'KnoxThom_Plot_' num2str(KT_paths) ...
    %     '_lay1_' num2str(lay1) '_lay2_' num2str(lay2) ...
    %     '_RA' num2str(RA) '_fromAsc_' num2str(rangeAscans(1)) '.png'], ...
    %     'Resolution',300);
    % **********************************************************************************

    % ----------- PLOT PHASE --------------------

    % ----------- PARAMETERS -----------
    N = Bscans;
    xx = xx(:);
    
    % ----------- COMPLEX CORRELATION MATRIX -----------
    P = complex(zeros(N,N));
    W = zeros(N,N);
    
    [d,i] = find(isfinite(m_phase) & isfinite(m_amp) & m_amp>0);
    j = i+d-1;
    valid = d>=2 & i<=N & j<=N;
    d = d(valid); i = i(valid); j = j(valid);
    
    idx = sub2ind(size(m_phase),d,i);
    z = exp(1i*m_phase(idx));
    w = m_amp(idx);
    
    P(sub2ind([N N],i,j)) = z;
    P(sub2ind([N N],j,i)) = conj(z);
    W(sub2ind([N N],i,j)) = w;
    W(sub2ind([N N],j,i)) = w;
    
    W = W./max(W(:));
    
    % ----------- GLOBAL PHASE SYNCHRONIZATION -----------
    H = W.*P;
    D = sum(W,2);
    A = spdiags(1./max(D,eps),0,N,N)*conj(H);
    
    [u,~] = eigs(A,1,'largestabs');
    u = u./max(abs(u),eps);
    u = u.*conj(u(1));
    
    % ----------- LOCAL PHASE INCREMENTS -----------
    Z = P(:,2:N).*conj(P(:,1:N-1));
    V = min(W(:,2:N),W(:,1:N-1));
    V = V.*(W(:,2:N)>0 & W(:,1:N-1)>0);
    
    C = sum(V.*Z,1);
    delta = angle(C(:));
    

    deltaGlobal = angle(u(2:end).*conj(u(1:end-1)));
    Rlocal = abs(C(:))./max(sum(V,1).',eps);
    
    % Detect ambiguous global phase transitions
    margin = 0.35;
    Rmin = 0.7;
    ambiguous = abs(deltaGlobal)>pi-margin;
    reliable = Rlocal>=Rmin & abs(delta)<pi/2;
    opposite = deltaGlobal.*delta<0;
    correct = ambiguous & reliable & opposite;
    
    % Select the branch using the local direction
    deltaFinal = deltaGlobal;
    deltaFinal(correct) = deltaGlobal(correct)-2*pi*sign(deltaGlobal(correct));
    
    % Reconstruct continuous phase
    colap1 = [0; cumsum(deltaFinal)];
    
    fprintf('Corrected transitions: %d\n',nnz(correct));
    
    % ----------- PHASE AGREEMENT -----------
    predicted = conj(P).*u.';
    error = angle(predicted.*conj(u));
    R = abs(sum(W.*exp(1i*error),2))./max(sum(W,2),eps);
    
    % ----------- PLOT -----------
    f2 = figure('Color','w','Position',[80 80 1200 650]);
    plot(xx,m_phase(:,column),'b-',xx,colap1,'r-','LineWidth',1.8);
    legend('Direct calculation (reference t_0 = 0)','Global phase with local branch selection','Location','northoutside','Orientation','horizontal','Box','off','FontName','Times New Roman','FontSize',15);
    set(gca,'FontName','Times New Roman','FontSize',18,'LineWidth',1.2,'Box','off','TickDir','out');
    grid on; xlim([xx(1) xx(end)]);
    xlabel('Time, t (s)','FontName','Times New Roman','FontSize',20);
    ylabel('Phase difference, \Delta\phi (rad)','FontName','Times New Roman','FontSize',20);
    

end