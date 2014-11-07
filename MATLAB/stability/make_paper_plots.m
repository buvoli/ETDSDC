%% === Set Parameters ===
clear; clc; close all;

% === Include common folder ===
wd = pwd();
cd ../common
addpath(pwd());
cd(wd);
% =============================

% Global Plot Parameters
plot_options      = struct( ...
    'plot_resolution' , 250, ...
    'lightest_contour', [.8,.8,.8], ...
    'darkest_contour' , [.5,.5,.5], ...
    'line_style'      , '-', ...
    'num_ticks'       , 7 , ...
    'print_title'     , false , ...
    'show_legend'      , false, ...
    'box'             , 'off');
label_colors      = {[.1 .1 .1], [53 146 177]/255, [253,209,7]/255, [214 7 37]/255};



% Global Save Parameters
save_options = struct( ...
    'renderer', 'zbuffer', ...
    'units'   , 'centimeters', ...
    'position', [0 0 8 7], ... %[x_position y_position width height]
    'dpi'     , 500);
plot_options.num_ticks = 5;


%% === 8th order ETD Methods ===

N = 8; M = N-1;
method = struct('type','ETD','N',N,'M',M);
base_filename = [method.type,'_',num2str(N), '_', num2str(M)];

% Stability Regions
plot_options.type     = 'stability';
plot_options.axis     = 40*[-1 1 -1 1];

filename = [base_filename,'_dissipative'];
r = -1*linspace(0,30,15);
plot_options.labels   = {{0,label_colors{2}}, {-15,label_colors{3}},{-30,label_colors{4}}}; % {r, color}
makePlot(method,r,plot_options);
set(gca, 'box','off');
savePlot(['stability/',filename],save_options);

filename = [base_filename,'_dispersive'];
%r = 1i*[0:2.5:25 25.5:.5:30];
r = 1i*linspace(0,30,15);
plot_options.labels   = {{0,label_colors{2}}, {15i,label_colors{3}},{30i,label_colors{4}}}; % {r, color}
makePlot(method,r,plot_options);
savePlot(['stability/',filename],save_options);

filename = [base_filename,'_both'];
r = exp(3*1i*pi/4)*linspace(0,30,15);
plot_options.labels   = {{0,label_colors{2}}, {15*exp(3*1i*pi/4),label_colors{3}},{30*exp(3*1i*pi/4),label_colors{4}}}; % {r, color}
makePlot(method,r,plot_options);
savePlot(['stability/',filename],save_options);

% Accuracy Regions
plot_options.type     = 'accuracy';
plot_options.axis     = 2*[-1 1 -1 1];
plot_options.epsilon  = 1e-8;

filename = [base_filename,'_dissipative'];
r = -1*linspace(0,5,15);
plot_options.labels   = {{0,label_colors{2}}, {-2.5,label_colors{3}},{-5,label_colors{4}}}; % {r, color}
makePlot(method,r,plot_options);
savePlot(['accuracy/',filename],save_options);

filename = [base_filename,'_dispersive'];
%r = 1i*([0:.25:2.25 2.5:.125:6]);
r = 1i*linspace(0,5,15);
plot_options.labels   = {{0,label_colors{2}}, {2.5i,label_colors{3}},{5i,label_colors{4}}}; % {r, color}
makePlot(method,r,plot_options);
savePlot(['accuracy/',filename],save_options);

filename = [base_filename,'_both'];
r = exp(3*1i*pi/4)*linspace(0,5,15);
plot_options.labels   = {{0,label_colors{2}}, {2.5*exp(3*1i*pi/4),label_colors{3}},{5*exp(3*1i*pi/4),label_colors{4}}}; % {r, color}
makePlot(method,r,plot_options);
savePlot(['accuracy/',filename],save_options);


%% === 8th order IMEX Methods ===

N = 8; M = N-1;
method = struct('type','IMEX','N',N,'M',M);
base_filename = [method.type,'_',num2str(N), '_', num2str(M)];

% Stability Regions
plot_options.type     = 'stability';
plot_options.axis     = 40*[-1 1 -1 1];

filename = [base_filename,'_dissipative'];
r = -1*linspace(0,30,15);
plot_options.labels   = {{0,label_colors{2}}, {-15,label_colors{3}},{-30,label_colors{4}}}; % {r, color}
makePlot(method,r,plot_options);
savePlot(['stability/',filename],save_options);

filename = [base_filename,'_dispersive'];
r = (1i)*linspace(0,30,15);
plot_options.labels   = {{0,label_colors{2}}, {15i,label_colors{3}},{30i,label_colors{4}}}; % {r, color}
makePlot(method,r,plot_options);
savePlot(['stability/',filename],save_options);

filename = [base_filename,'_both'];
r = exp(3*1i*pi/4)*linspace(0,30,15);
plot_options.labels   = {{0,label_colors{2}}, {15*exp(3*1i*pi/4),label_colors{3},'5e^{3i\pi/4}'},{30*exp(3*1i*pi/4),label_colors{4},'10e^{3i\pi/4}'}}; % {r, color, label?}
makePlot(method,r,plot_options);
savePlot(['stability/',filename],save_options);

% Accuracy Regions
plot_options.type     = 'accuracy';
plot_options.axis     = 2*[-1 1 -1 1];
plot_options.epsilon  = 1e-8;

filename = [base_filename,'_dissipative'];
r = -1*linspace(0,3,22);
plot_options.labels   = {{0,label_colors{2}}, {-1,label_colors{3}},{-2,label_colors{4}}}; % {r, color}
makePlot(method,r,plot_options);
savePlot(['accuracy/',filename],save_options);

filename = [base_filename,'_dispersive'];
r = 1i*linspace(0,3,22);
plot_options.labels   = {{0,label_colors{2}}, {1i,label_colors{3}},{2i,label_colors{4}}}; % {r, color}
makePlot(method,r,plot_options);
savePlot(['accuracy/',filename],save_options);

filename = [base_filename,'_both'];
r = exp(3*1i*pi/4)*linspace(0,3,22);
plot_options.labels   = {{0,label_colors{2}}, {exp(3*1i*pi/4),label_colors{3}},{2*exp(3*1i*pi/4),label_colors{4}}}; % {r, color}
makePlot(method,r,plot_options);
savePlot(['accuracy/',filename],save_options);

%% === 16th order ETD Methods ===

N = 16; M = N-1;
method = struct('type','ETD','N',N,'M',M);
base_filename = [method.type,'_',num2str(N), '_', num2str(M)];

% Stability Regions
plot_options.type     = 'stability';
plot_options.axis     = 40*[-1 1 -1 1];

filename = [base_filename,'_dissipative'];
r = -1*linspace(0,30,15);
plot_options.labels   = {{0,label_colors{2}}, {-15,label_colors{3}},{-30,label_colors{4}}}; % {r, color}
makePlot(method,r,plot_options);
savePlot(['stability/',filename],save_options);

filename = [base_filename,'_dispersive'];
r = 1i*[linspace(0,30,15) 60];
plot_options.labels   = {{0,label_colors{2}}, {15i,label_colors{3}},{30i,label_colors{4}}, {60i, label_colors{1}}}; % {r, color}
makePlot(method,r,plot_options);
savePlot(['stability/',filename],save_options);

filename = [base_filename,'_both'];
r = exp(3*1i*pi/4)*linspace(0,30,15);
plot_options.labels   = {{0,label_colors{2}}, {15*exp(3*1i*pi/4),label_colors{3}},{30*exp(3*1i*pi/4),label_colors{4}}}; % {r, color}
makePlot(method,r,plot_options);
savePlot(['stability/',filename],save_options);

% Accuracy Regions
plot_options.type     = 'accuracy';
plot_options.axis     = 16*[-1 1 -1 1];
plot_options.epsilon  = 1e-8;

filename = [base_filename,'_dissipative'];
r = -1*linspace(0,30,15);
plot_options.labels   = {{0,label_colors{2}}, {-15,label_colors{3}},{-30,label_colors{4}}}; % {r, color}
makePlot(method,r,plot_options);
savePlot(['accuracy/',filename],save_options);

filename = [base_filename,'_dispersive'];
r = 1i*[linspace(0,15,15) 16:18];
plot_options.labels   = {{0,label_colors{2}}, {7.5i,label_colors{3}},{15i,label_colors{4}}}; % {r, color}
makePlot(method,r,plot_options);
savePlot(['accuracy/',filename],save_options);

filename = [base_filename,'_both'];
r = exp(3*1i*pi/4)*linspace(0,30,15); r(7) = 13*exp(3*1i*pi/4); r(13) = 26*exp(3*1i*pi/4);
plot_options.labels   = {{0,label_colors{2}}, {13*exp(3*1i*pi/4),label_colors{3}},{26*exp(3*1i*pi/4),label_colors{4}}}; % {r, color}
makePlot(method,r,plot_options);
savePlot(['accuracy/',filename],save_options);

%% === 16th order IMEX Methods ===

N = 16; M = N-1;
method = struct('type','IMEX','N',N,'M',M);
base_filename = [method.type,'_',num2str(N), '_', num2str(M)];

% Stability Regions
plot_options.type     = 'stability';
plot_options.axis     = 40*[-1 1 -1 1];

filename = [base_filename,'_dissipative'];
r = -1*linspace(0,30,15);
plot_options.labels   = {{0,label_colors{2}}, {-15,label_colors{3}},{-30,label_colors{4}}}; % {r, color}
makePlot(method,r,plot_options);
savePlot(['stability/',filename],save_options);

filename = [base_filename,'_dispersive'];
r = (1i)*linspace(0,30,15);
plot_options.labels   = {{0,label_colors{2}}, {15i,label_colors{3}},{30i,label_colors{4}}}; % {r, color}
makePlot(method,r,plot_options);
savePlot(['stability/',filename],save_options);

filename = [base_filename,'_both'];
r = exp(3*1i*pi/4)*linspace(0,30,15);
plot_options.labels   = {{0,label_colors{2}}, {15*exp(3*1i*pi/4),label_colors{3}},{30*exp(3*1i*pi/4),label_colors{4}}}; % {r, color}
makePlot(method,r,plot_options);
savePlot(['stability/',filename],save_options);

% Accuracy Regions
plot_options.type     = 'accuracy';
plot_options.axis     = 16*[-1 1 -1 1];
plot_options.epsilon  = 1e-8;

filename = [base_filename,'_dissipative'];
r = -1*linspace(0,20,15);
plot_options.labels   = {{0,label_colors{2}}, {-10,label_colors{3}},{-20,label_colors{4}}}; % {r, color}
makePlot(method,r,plot_options);
savePlot(['accuracy/',filename],save_options);

filename = [base_filename,'_dispersive'];
r = 1i*linspace(0,9,15);
plot_options.labels   = {{0,label_colors{2}}, {4.5i,label_colors{3}},{9i,label_colors{4}}}; % {r, color}
makePlot(method,r,plot_options);
savePlot(['accuracy/',filename],save_options);

filename = [base_filename,'_both'];
r = exp(3*1i*pi/4)*[linspace(0,16,15) 16+16/14:16/14:20];
plot_options.labels   = {{0,label_colors{2}}, {8*exp(3*1i*pi/4),label_colors{3}},{16*exp(3*1i*pi/4),label_colors{4}}}; % {r, color}
makePlot(method,r,plot_options);
savePlot(['accuracy/',filename],save_options);
