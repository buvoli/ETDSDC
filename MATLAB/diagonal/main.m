% =====================================================================
% main.m runs each numerical method on a given equation with specificed
% timestep and plots results for each method, as well as relative error
% with respect to the ETDSDC method.
% =====================================================================

clear; clc; close all;

% === Include common folder ===
wd = pwd();
cd ../common
addpath(pwd());
cd(wd);
% =============================

% Load Equation
equation_name = 'kuramoto';  % kuramoto, kdv, nls, nikolaevskiy, quasigeostrophic

    cd(['equations/',equation_name]);
    init
    cd ../../

% Set options for each method
rk_options      = struct('parameters',pars,'max_ts_to_store',100);
etdsdc_options  = struct('n',8,'m',7,'parameters',pars,'max_ts_to_store',100);
imsdc_options   = struct('n',8,'m',4,'parameters',pars,'max_ts_to_store',100);

% Set number of timesteps
Nt_edc  = 300;
Nt_idc  = 2000;
Nt_rk   = 20000;

% Run Methods
[~, YEDC, tcpu_Edc, tccpu_Edc] = etdsdc(LF,NF,tspan,y0,Nt_edc,etdsdc_options); 
[~, YIDC, tcpu_Idc, tccpu_Idc] = imexsdc(LF,NF,tspan,y0,Nt_idc,imsdc_options);
[~, YRK,  tcpu_rk,  tccpu_rk]  = etdrk4(LF,NF,tspan,y0,Nt_rk,rk_options);

% Plot Results
scrsz = get(groot,'ScreenSize');
figure_width = 1000; figure_height = 400;
figure('Position',[((scrsz(3)-figure_width)/2) ((scrsz(4)-figure_height)/2) figure_width figure_height])

if(equation_dimension == 1)
    subplot(1,3,1); mesh(filter(YEDC)); axis tight; title('ETDSDC Solution');  set(gca,'Ydir','reverse'); xlabel('x'); ylabel('t'); view([-90 90]);
    subplot(1,3,2); mesh(filter(YIDC)); axis tight; title('IMEXSDC Solution'); set(gca,'Ydir','reverse'); xlabel('x'); ylabel('t'); view([-90 90]);
    subplot(1,3,3); mesh(filter(YRK));  axis tight; title('ETDRK4 Solution');  set(gca,'Ydir','reverse'); xlabel('x'); ylabel('t'); view([-90 90]);
elseif(equation_dimension == 2)
    subplot(1,3,1); mesh(ys,xs,filter(YEDC(:,end))); axis tight; title(['ETDSDC Solution at t=',num2str(tspan(end))]);  xlabel('x'); ylabel('y'); view([-90 90]);
    subplot(1,3,2); mesh(ys,xs,filter(YIDC(:,end))); axis tight; title(['ETDSDC Solution at t=',num2str(tspan(end))]);  xlabel('x'); ylabel('y'); view([-90 90]);
    subplot(1,3,3); mesh(ys,xs,filter(YRK(:,end)));  axis tight; title(['ETDSDC Solution at t=',num2str(tspan(end))]);  xlabel('x'); ylabel('y'); view([-90 90]);
end
colormap bone;

% Errors (Relative to ETDSDC)
error_1 = error_filter(YRK(:,end),YEDC(:,end));  m1 = mean(error_1);
error_2 = error_filter(YIDC(:,end),YEDC(:,end)); m2 = mean(error_2);

fprintf('EDC\t\tTCPU: %.4f\t TCCPU: %.4f\t ME: %.4e\n', tcpu_Edc, tccpu_Edc, 0);
fprintf('RK4\t\tTCPU: %.4f\t TCCPU: %.4f\t ME: %.4e\n', tcpu_rk, tccpu_rk, m1);
fprintf('IDC\t\tTCPU: %.4f\t TCCPU: %.4f\t ME: %.4e\n', tcpu_Idc, tccpu_Idc, m2);