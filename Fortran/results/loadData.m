% load_data imports results from a Fortran Run
clear; close all; clc;

% kuramoto, nikolaevskiy, quasigeostrophic, kdv, nls
eqn_name = 'nls';
save_dir = 'figures/';

save_options = struct(              ...
    'renderer', 'painters',         ...
    'units'   , 'centimeters',      ...
    'position', [0 0 40 22],        ... %[x_position y_position width height]
    'dpi'     , 200);

plot_options      = struct(         ...
    'MarkerSize' , 4,               ...
    'lineWidth'  , 1,               ...
    'titleFontWeight', 'bold',      ...
    'fontSize', 11,                 ...
    'orderLines', true);
line_colors = struct('o4',[53 146 177]/255,'o8',[253,209,7]/255,'o16',[214 7 37]/255,'o32',[78 128 72]/255); %'o7',[153,143,184]/255
line_styles = struct('imexsdc','o-','etdsdc','d-','etdrk4','s-');

% Load data arrays
Nts     = load([eqn_name,'/Nts.txt']);
hs      = load([eqn_name,'/hs.txt']);
times   = load([eqn_name,'/times.txt']);
error   = load([eqn_name,'/errors.txt']);
Fs      = round(load([eqn_name,'/Fs.txt']));

% === Note: To load complex data ===
%RD   = load([eqn_name,'/complex_data.txt']);
%data = RD(:,1:2:end) + 1i*RD(:,2:2:end);

line_colors = struct('o4',[53 146 177]/255,'o8',[253,209,7]/255,'o16',[214 7 37]/255,'o32',[78 128 72]/255);
line_styles = struct('imexsdc','o-','etdsdc','d-','etdrk4','s-');

% Equation parameters
etd_methods  = [4 8 16 32]; num_etd = length(etd_methods);
imex_methods = [4 8 16 32]; num_imex = length(imex_methods);
tspan = [0 30];
equation_dimension = 1;

% === Generate Visual Solution (using MATLAB code) ===
matlab_implementation_path = '../../MATLAB';
eqn_type = 'diagonal';
matlab_eqn_dir = [matlab_implementation_path, '/', eqn_type,'/equations/',eqn_name];
cur_dir    = pwd();

if(exist(matlab_eqn_dir,'dir'))
    
    % === Include common folder ===
    wd = pwd();
    cd([matlab_implementation_path, '/common']);
    addpath(pwd());
    cd(wd);
    % =============================
   
    % Solve Equation
    cd(matlab_eqn_dir)
    init
    cd ../../
    etdsdc_options = struct('n',4,'m',3,'parameters',pars,'max_ts_to_store',300);
    min_ind = find(error(:,2)<1e-3,1);
    Nt = max(ceil(Fs(min_ind)/(etdsdc_options.n * (etdsdc_options.m+1))),300);
    [T, Y] = etdsdc(LF,NF,tspan,y0,Nt,etdsdc_options);
    cd(cur_dir)
else
    filter = @(X) X;
    xs = 0; ys = 0; Y = 0; T = 0; solution_name = 'na';
    fprintf(['WARNING: Solution Figure was not generated. \n directory not found - ', matlab_eqn_dir,'\n']); 
end


% === Create Plots ====
figure(1); set(gcf, 'Position', [0 0 1000 700]);
% subplot arrangment (2x2 grid + legend on right)
right_border    = 0;
subplot_margin  = .1;
subplot_width   = (1-right_border - 3*subplot_margin)/2;
subplot_height  = (1-3*subplot_margin)/2;

s1 = subplot('Position',[.1 .2+subplot_height subplot_width subplot_height]);
s2 = subplot('Position',[.2+subplot_width .2+subplot_height subplot_width subplot_height]);
s3 = subplot('Position',[.1 .1 subplot_width subplot_height]); 
s4 = subplot('Position',[.2+subplot_width .1 subplot_width subplot_height]);

% Plot Visual Solution
if(equation_dimension == 1)
    axes(s1); 
    %surf(T,xs,filter(Y),'EdgeColor', 'none'); view([-90 90]); colormap bone; set(s1,'ydir','reverse','xlim',[tspan(1) tspan(end)],'ylim',[xs(1) xs(end)]);
    imagesc(xs,T,transpose(filter(Y))); set(gca,'YDir','normal'); colormap bone;
    title('Solution Plot','fontweight',plot_options.titleFontWeight,'fontsize',plot_options.fontSize); ylabel('t'); xlabel('x');
elseif(equation_dimension == 2)
    axes(s1); 
    %surf(xs,ys,filter(Y(:,end)),'EdgeColor', 'none'); view([0 90]); colormap bone; set(s1,'ylim',[ys(1) ys(end)],'xlim',[xs(1) xs(end)]);
    imagesc(ys,xs,filter(Y(:,end))); set(gca,'YDir','normal'); colormap bone;
    title([solution_name,' Plot at Time t=',num2str(tspan(end))],'fontweight','bold','fontsize',12); xlabel('x'); ylabel('y');
end    

% Build Legend
dummy_x = [0 1]; dummy_y = [1e-100,1e-100]; 
legend_names = {}; fields = fieldnames(line_colors); axes(s2);
for j=1:length(fields)
    loglog(dummy_x, dummy_y, '-', 'lineWidth',plot_options.lineWidth,'color',line_colors.(fields{j}));  hold on;
    legend_names{end+1} = [fields{j}(2:end), 'th order']; 
end
fields = fieldnames(line_styles);
for j=1:length(fields)
    loglog(dummy_x,dummy_y,line_styles.(fields{j}),'MarkerFaceColor',[1 1 1], 'MarkerEdgeColor',[0 0 0], 'MarkerSize',plot_options.MarkerSize, 'lineWidth', plot_options.lineWidth, 'color',[.99 .99 .99])
    legend_names{end+1} = fields{j};
end
hL = legend(legend_names);
set(hL,'Location', 'SouthWest','box','off');
% End Legend


% Order Lines
ol = zeros(size(hs));
for j=1:(num_etd + num_imex + 1)
    if(j <= num_etd)
        order = (etd_methods(j));
    elseif(j <= (num_etd + num_imex))
        order = (imex_methods(j - num_etd));
        break;
    else
        order = 4;
        break;
    end
    
    % Fit order line
    error_tol = 10*min(error(:)); % only plot lines above 1e-10
    nne = find(~isnan(error(:,j)) == 1,1);
    if(error(nne,j) > error_tol)
        cc  = error(nne,j)/(hs(nne,j)^order);
        ol(:,j) = [NaN*ones(nne-1,1); cc*hs(nne:end,j).^order];
    end
end
if(plot_options.orderLines)
    axes(s3); loglog(hs,ol,'--','color',[.5 .5 .5]); hold on;
end

% Plot Methods
for j=1:(num_etd + num_imex + 1)

    if(j <= num_etd)
        order = num2str(etd_methods(j));
        name  = 'etdsdc';
    elseif(j <= (num_etd + num_imex))
        order = num2str(imex_methods(j - num_etd));
        name  = 'imexsdc';
    else
        order = '4';
        name  = 'etdrk4';
    end

    axes(s2); loglog(Fs, error(:,j),line_styles.(name),'MarkerEdgeColor',line_colors.(['o',order])/1.5,'MarkerFaceColor',[1 1 1], 'MarkerSize',plot_options.MarkerSize, 'lineWidth', plot_options.lineWidth,'color',line_colors.(['o',order]));  hold on;
        title('Accuracy vs Function Evaluations','fontweight',plot_options.titleFontWeight,'fontsize',plot_options.fontSize); 
        xlabel('Function Evaluations','fontsize',plot_options.fontSize);
        ylabel(['Relative Error at t=',num2str(tspan(end))],'fontsize',plot_options.fontSize)
     axes(s3); loglog(hs(:,j),error(:,j), line_styles.(name),'MarkerEdgeColor',line_colors.(['o',order])/1.5,'MarkerFaceColor',[1 1 1], 'MarkerSize',plot_options.MarkerSize, 'lineWidth', plot_options.lineWidth,'color',line_colors.(['o',order])); hold on;
        title('Accuracy vs Stepsize','fontweight',plot_options.titleFontWeight,'fontsize',plot_options.fontSize);
        xlabel('Step Size h','fontsize',plot_options.fontSize);
        ylabel(['Relative Error at t=',num2str(tspan(end))],'fontsize',plot_options.fontSize);
     axes(s4); loglog(times(:,j), error(:,j),line_styles.(name),'MarkerEdgeColor',line_colors.(['o',order])/1.5,'MarkerFaceColor',[1 1 1], 'MarkerSize',plot_options.MarkerSize, 'lineWidth', plot_options.lineWidth,'color',line_colors.(['o',order])); hold on;
        title('Accuracy vs Computer Time','fontweight',plot_options.titleFontWeight,'fontsize',plot_options.fontSize);
        xlabel('Computer Time (sec)','fontsize',plot_options.fontSize);
        ylabel(['Relative Error at t=',num2str(tspan(end))],'fontsize',plot_options.fontSize);
end

sf = 1.75;
NNE = ~isnan(error);
max_nonnan_h    = max(hs(NNE));   min_nonnan_h    = min(hs(NNE));   
max_nonnan_time = max(times(NNE)); min_nonnan_time = min(times(NNE));
for i=1:(num_etd + num_imex)
    COLS = or(NNE(:,i),NNE(:,i+1));
end
max_nonan_F = max(Fs(COLS)); min_nonan_F = min(Fs(COLS));

set(s1,'box','off')
set(s2,'xlim',[(1/sf)*min(Fs(:)) (sf)*max(Fs(:))], 'ylim', [1e-12 1e1], 'YTick',10.^(-(12:-2:0)), 'box', 'off')
set(s3,'xlim',[(1/sf)*min_nonnan_h (sf)*max_nonnan_h], 'ylim', [1e-12 1e1], 'YTick',10.^(-(12:-2:0)), 'box', 'off')
set(s4,'xlim',[(1/sf)*min_nonnan_time (sf)*max_nonnan_time], 'ylim', [1e-12 1e1], 'YTick',10.^(-(12:-2:0)), 'box', 'off') 

savePlot([save_dir,eqn_name], save_options);