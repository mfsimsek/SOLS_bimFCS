clear;
warning('off', 'all');
%AC Colors LUT
colors=[0.09,0.96,0.04;0.26	0.54,0.88;0.95,0.05,0.05;0.03,0.23,0.79;...
0.92,0.48,0.09;0.73	0.62,0.26;0.48,0.67,0.33;0.57,0.39,0.67;...
0.23,0.36,0.13;0.45	0.98,0.72;0.19,0.96,0.04;0.36,0.54,0.88;...
0.05,0.05,0.05;0.13	0.23,0.79;0.02,0.48,0.09;0.83,0.62,0.26;...
0.58,0.67,0.33;0.67	0.39,0.67;0.33,0.36,0.13;0.55,0.98,0.72];
% Update input values in between ----------------------------
FileTif = 'sample_e3.ome.tif'; %name of the tif file to analyze

frame_time=0.008;% frame duration in sec. Confirm from metadata.
pixel_size=0.39; % camera pixel size in um
sigma= 387.95; %PSF waist in nm. Measured using sub-diffraction limit sized beads.
bg= 104; %camera dark noise average/pixel

minbin=2; % Minimum super-pixel to anaylze. If signal is low, begin with bin 2 or 3. Bin1 data will be noisy.
maxbin=8;  % Maximum super-pixel to anaylze. Check considerations in Weixiang, Simsek, Pralle, Methods, 2018. 8-10 for bin2, 12-16 for bin1 acquisition is a good start.

first_frame= 2000; % skip first t
last_frame=20000; % exclude last t (or time-split analysis)

%for z= 0:3 % for space-split analysis %activate last line as well
%spacer=100; %space-split analysis width
firstx= 0; % skip leftmost x pixels
firsty= 10; % skip topmost y pixels (3-10 for camera artifacts)
%firsty= 10+z*spacer/2; % skip topmost y pixels + space-split analysis
lastx= 110; % exclude right (x) pixels
lasty= 260; % exclude bottom (y) pixels
%lasty= spacer+firsty+z*spacer/2; % space-split analysis

% --- w2_bin for light sheet channels bin2 acquisition ---
%515nm Laser channel
w2_bin=[0.661417755733266 0.830532703104362	1.14395226374891	1.61633098064579	2.23880573677181	2.99055449341354	3.85420229112375	4.82376960607463	5.89407341217341	7.06502016851100	8.33712616923030	9.70907384608769	11.1809734825458	12.7546486390940	14.4278258228198]; %effective beam waists for square convoluted profile as calculated at beamwaister_FCS.m
% Update input values in between ----------------------------

%=====Mise-en-scene================
w2_bin(maxbin+1:end)=[];
w2_bin(1:minbin-1)=[];
bins=maxbin-minbin+1;
%==================================
%read image file into matrix
%==================================
InfoImage = imfinfo(FileTif);
s=InfoImage.Filename;
s1=strsplit(s,'\');
s1=strsplit(s1{size(s1,2)},'.');
s=s1{1};

% choose class based on BitDepth
if InfoImage(1).BitDepth == 8
    cls = 'uint8';
elseif InfoImage(1).BitDepth == 16
    cls = 'uint16';
else
    cls = 'double';
end
% Preallocate a 4D matrix for speed
Image = zeros(InfoImage(1).Height, InfoImage(1).Width, length(InfoImage), cls); 
TifLink = Tiff(FileTif, 'r');
for i = 1:length(InfoImage)
    TifLink.setDirectory(i);
    Image(:,:,i) = TifLink.read();
end
TifLink.close();

stack = Image(firsty+1:lasty,firstx+1:lastx,first_frame+1:last_frame);

[h, w, numFrames] = size(stack);
    
    fprintf('\n=== Light-Sheet FCS Analysis with Diffusion Model Fitting ===\n');
    fprintf('Stack size: %dx%d x %d frames\n', h, w, numFrames);
    fprintf('Frame time: %.4f s | Pixel size: %.2f um\n\n', frame_time, pixel_size);

s=strcat(s,"fi_",num2str(first_frame/1000),'_',num2str(numFrames/1000),"K_x",num2str(firstx),"_",num2str(lastx),"_y",num2str(firsty),"_",num2str(lasty));
diary(strcat('FCS_outputlog_',s,'.txt')); %create diary of analysis

%==================================
%bleach correction
%==================================

fprintf('Computing bleach correction...\n');
    
    stack_dbl = double(stack);
    time_vector = (0:numFrames-1)' * frame_time;
    
    % global trajectory for spatial average of pixels
    stack_mean = mean(stack_dbl, [1, 2]);
    stack_means=stack_mean.*ones(h,w,numFrames);
    stack_mean = squeeze(stack_mean);
    
    stack_corrected = stack_dbl-stack_means+stack_mean(1);

%==================================
%lin to log resampling / autocorrelate
%================================== 
fprintf('Computing auto-correlations...\n');
max_lag_idx = floor(numFrames / 2);
tau_lin = (1:max_lag_idx) * frame_time;

AC=cell(1,bins);
lags=cell(1,bins);
for b=minbin:maxbin
f_bin=bin2d_overlap(stack_corrected,b)- bg; %creates super-pixelated data with 2d span array in bin2d_overlap.m
h=size(f_bin,1);
w=size(f_bin,2);
G_lin=zeros(h,w, (max_lag_idx*2)+1);
f_bincorr=f_bin;
for i=1:h
    for j=1:w
        f_bincorr(i,j,:)=f_bin(i,j,:)/mean(f_bin(i,j,:))-1;
        G_lin(i,j,:)=xcorr(f_bincorr(i,j,:),max_lag_idx);
    end
end
[AC{b-minbin+1},lags{b-minbin+1}]=lintolog(mean(mean(G_lin))); %uses lintolog.m file
end

%==================================
%FCS function double comp fitting / linearized Nfast
%==================================
fprintf('Fitting log-resampled auto-correlation functions...\n');
clear lagtime2c AC2c
%bins=3; use this individually if linearizing Nfast over smallest 3 super-pixels
for j=1:bins
    lagtime2c{j} = lags{j}(~isnan(lags{j}))*frame_time;
    AC2c{j} = AC{j}(~isnan(AC{j}));
end
AC2c=cellfun(@(x) x(1:85), AC2c, 'UniformOutput',false);
lagtime2c=cellfun(@(x) x(1:85), lagtime2c, 'UniformOutput',false);
    

[param,y_Fit]=createFit_2c(lagtime2c,AC2c, w2_bin, sigma, pixel_size, minbin);
yFit=cell2mat(y_Fit)';
lagtime=lagtime2c{1};

%==================================
%plot G(tau) vs. tau with fits
%==================================
Gvstau=figure;
semilogx(lagtime, AC2c{1},'Color',colors(1+minbin-1,:),'Marker','o');
hold on;
semilogx(lagtime, squeeze(yFit(1,:)),'Color',colors(1+minbin-1,:),'LineStyle','-');
for p=2:bins
semilogx(lagtime, AC2c{p}, 'Color',colors(p+minbin-1,:),'Marker','o');
semilogx(lagtime, squeeze(yFit(p,:)),'Color',colors(p+minbin-1,:),'LineStyle','-'); % Plot data as markers, fit as a line
end
hold off;

% Label axes
xlabel( 'lagtime(s)');
ylabel( 'G(tau)');

set(Gvstau, 'color', 'none'); 
set(gcf, 'color', 'none'); 
set(gca, 'color', 'none');
set(gca, "YLim", [0,prctile(AC{1},98)*1.2]);
% Export the figure with a transparent background
exportgraphics(Gvstau, strcat('2C_ACplot_bin_',s,'.png'), 'BackgroundColor', 'none');

%==================================
%FCS law plot
%==================================
fprintf('Performing FCS law calculations...\n');
% Preallocate arrays to prevent them from "growing" during the loop
Tau_bin2 = zeros(1, bins);
%Tau_bin2_err = zeros(1, bins);
C_bin2 = zeros(1, bins);
%C_bin2_err = zeros(1, bins);
Cn_bin2 = zeros(1, bins);
%Cn_bin2_err = zeros(1, bins);
Area = zeros (1, bins);
RatioN = zeros (1, bins);

for b=1:bins
%    cfV = coeffvalues(fitresults{b});
    Df=param(2)*1e-6;   %----- for nm^2 to um^2
    Nf = w2_bin(b)*1e6/(param(1));  %----- fast comp
    Ns = 1/param((b-1)*3 + 4); %----- slow comp
    Ds = param((b-1)*3 + 3)*1e-6;   %----- for nm^2 to um^2
    Dapp=(Nf*Df+Ns*Ds)/(Nf+Ns);  %----- apparent diffusion coeff
    N=Nf+Ns; %----- total diffusers
    RatioN(b)=Nf/N; %----- ratio of fast to total population

    Tau_bin2(b)=1000*(w2_bin(b)/(4*Dapp));
    Area(b) = (sigma*(exp(-(b+minbin-1)^2*pixel_size^2*250000/sigma^2)-1)/500/sqrt(pi)/(b+minbin-1)^2/pixel_size^2+erf((b+minbin-1)*pixel_size*500/sigma)/(b+minbin-1)/pixel_size)^-2;
    C_bin2(b) = 10*N*numFrames/(Area(b)*140*6.02); % in nM, calibration constant comes from known concentration bead data.
    Cn_bin2(b) = (4*pi()*(2.8)^3/3)*0.602*C_bin2(b); % per cell count for 5.6 micron internuclear distance as measured in PSM tissue
end

% Prepare FCS law data for the linear fit function.
x_data = w2_bin(:);
y_data = Tau_bin2(:);

j0=0;
for j=1:20
    if (prctile(Tau_bin2,5*j,2)>2*prctile(Tau_bin2,5*j-5,2))
        rend=2*prctile(Tau_bin2,5*j-5,2);
        j0=1;
    end
end
if (j0==0)
    rend=max(Tau_bin2);
end
tf=excludedata(w2_bin(:),Tau_bin2(:),'Range',[0 rend]);

fcs=figure; % Create a new figure window
e = plot(w2_bin(~tf), Tau_bin2(~tf), 'o', 'DisplayName', 'Data');
hold on; % Keep the current plot active to add more elements

%==================================
% Calculate and Plot the Linear Fit
%==================================

% Generate y-values for the fitted line over a smooth range
x_fit = linspace(0, max(w2_bin), 100)'; % Smooth range for the fit line

mdl=fitlm(w2_bin(~tf), Tau_bin2(~tf), 'Linear', 'RobustOpts','on'); %linear model creation
y_fit=predict(mdl,x_fit); %fit from the linear model
% Overlay the fitted line
plot(x_fit, y_fit, '-r', 'LineWidth', 1.5, 'DisplayName', 'Linear Fit'); % Use '-r' for a red line

% Label axes
xlabel( 'Area (\mum^2)');
ylabel( 'Transit time (ms)');
grid on

set(fcs, 'color', 'none'); 
set(gcf, 'color', 'none'); 
set(gca, 'color', 'none');
ylim ([0 median(Tau_bin2)*3]); %limits of y-axis in display

% Export the figure with a transparent background
exportgraphics(fcs, strcat('2C_FCS_law_',s,'.png'), 'BackgroundColor', 'none');

b=table2array(mdl.Coefficients(1,1));
m=table2array(mdl.Coefficients(2,1));
D_eff=250/m; %effective diffusion coefficient from FCS law fit

fcslawstr=sprintf('t_0: %.2f ms | m: %.2f | D_eff: %.2f um^2/s\n\n', b, m, D_eff); %fit results string

%===========================
%Starting from bin5 (after initial non-linear trend)
%===========================
fcs2=figure;
plot(w2_bin,Tau_bin2,'bo');
hold on
mdl=fitlm(w2_bin(5:end), Tau_bin2(5:end), 'Linear', 'RobustOpts','on');
y_fit=predict(mdl,x_fit);
% Overlay the fitted line
plot(x_fit, y_fit, '-r', 'LineWidth', 1.5, 'DisplayName', 'Linear Fit'); % Use '-r' for a red line
hold off

exportgraphics(fcs2, strcat('2C_after5_FCS_law_',s,'.png'), 'BackgroundColor', 'none');
b=table2array(mdl.Coefficients(1,1));
m=table2array(mdl.Coefficients(2,1));
D_eff=250/m;

%Measure the confinement for 2-component difffusive behaviour
yp=1000;
xint=0;
while abs(predict(mdl,xint))>yp
    xint=xint+0.01;
    yp=predict(mdl,xint);
end
L = sqrt(xint*pi()/2);
Sconf=Tau_bin2.*4e-9*param(1)./(w2_bin-xint);

fcslawstr2=sprintf('t_0: %.2f ms | m: %.2f | D_eff: %.2f um^2/s | L: %.2f um \n\n', b, m, D_eff, L);

%==================================
% Create a table
%==================================

results=[C_bin2(:), Cn_bin2(:), RatioN(:), Tau_bin2(:), w2_bin(:), Sconf(:)];
col_headers = {'C(nM)', 'Cn(nuclei)', 'Ratio(f)', 'Tau(ms)', 'Area', 'Sconf'};
T = array2table(results, 'VariableNames', col_headers);
writetable(T, strcat('2CimFCSresults_',s,'.xlsx'));
T2 = array2table(param);
writetable(T2, strcat('2CimFCSresults_',s,'.xlsx'), 'Sheet', 2);
writematrix(fcslawstr, strcat('2CimFCSresults_',s,'.xlsx'), 'Sheet', 3, 'Range', 'A1');
writematrix(fcslawstr2, strcat('2CimFCSresults_',s,'.xlsx'), 'Sheet', 3, 'Range', 'A2');
diary off;
%end %activate this for space-split anaysis