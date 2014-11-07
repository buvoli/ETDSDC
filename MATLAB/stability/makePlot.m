function makePlot(method,rs,plot_options)
%MAKEPLOT Generates a figure with overalyed 2D stability regions for ETDSDC_N^M
%or IMEXSDC_N^M schemes assuming fixed r=h\nu and chebyshev quadrature points.
% Parameters
%   method      (struct) Contains method parameters. Fields:
%                 type (string)  : ETD or IMEX 
%                 N    (integer) : Number of quadrature points
%                 M    (integer) : Number of correction sweeps
%
%   rs          (array) values of r=h\nu to calculate stability region for.
%
%   plot_option (struct) contains plotting options. Fields:
%                 type                (string)  : 'accuracy' or 'stability'.
%                 epsilon             (double)  :  epsilon value for accuracy plot.
%                 labels              (cellarray) : Contains r values to be labeled. NOTE: labels{i}={r_value, color}.
%                 axis                (array)   : [min_x max_x min_y right_y].
%                 lightest_countour   (array)   : color for first r value.
%                 darkest_countour	  (array)   : color for largest r value.
%                 plot_resolution     (integer) :  number of points to plot in each dimension.
%                 line_style          (string)  : line style for contours.
%                 num_ticks           (integer) : number of ticks per axis.
%                 print_title         (boolean) : shows/hides title.
%                 show_legend         (boolean) : shows/hides legend.
%                 box                 (boolean) : shows/hides box around graph.
%                   
% Returns
%   N/A

% Plot Parameters
r_to_label = cellfun(@(x) x{1}, plot_options.labels);
delta      = (plot_options.lightest_contour - plot_options.darkest_contour)/(length(rs)-1);
handles    = cell(0);
handle_count = 1;
lx = abs(plot_options.axis(1));
ly = abs(plot_options.axis(3));
Nxy = plot_options.plot_resolution;



figure(1); clf;
for k=1:length(rs)
    r = rs(k);
    
    % set chebyshev quadrature points
    tau = 0.5 - 0.5*cos(pi*((0:method.N-1)')/(method.N-1));
    
    % Initialize Integration Matrix
    if(strcmp(method.type,'ETD'))
        amp = @ampETD;
        I = initW(r,1,tau);
    else
        amp = @ampIMEX;
        I = initIM(1,tau);  % divide by two since integration interval is [-1 1] instead of [0 1]
    end
    
    % Calculate Stability/Accuracy Region
    a = linspace(-lx,lx,Nxy);
    b = linspace(-ly,ly,Nxy);
    A = zeros(Nxy,Nxy);
    parfor i=1:Nxy
        for j=1:Nxy
            z       = (a(i)+b(j)*1i);
            psi = amp(method.N,method.M,r,z,I);            
            if(strcmp(plot_options.type,'accuracy'))
                A(i,j)  = abs(psi - exp(r+z));
            elseif(strcmp(plot_options.type,'stability'))
                A(i,j)  = abs(psi);
            end
        end
    end
    
    % Select Line color/size
    i = find(r_to_label == r,1);
    if(~isempty(i))
       line_color   = plot_options.labels{i}{2};
       line_size    = 1.3; 
       store_handle = true;
    else
       line_color   = plot_options.darkest_contour+(k-1)*delta;
       line_size    = 0.5;
       store_handle = false;
    end
    
    % Plot Contours
    if(strcmp(plot_options.type,'accuracy'))
        [~,h] = contour(b,a,A,[plot_options.epsilon plot_options.epsilon],'color',line_color,'LineWidth',line_size,'LineStyle',plot_options.line_style);
    elseif(strcmp(plot_options.type,'stability'))
       [~,h] = contour(b,a,A,[1 1],'color',line_color,'LineWidth',line_size,'LineStyle',plot_options.line_style);
    end
    
    if(store_handle)
        handles{handle_count} = h;
        handle_count = handle_count +1;
    end
    drawnow;
    hold on;
end

% place colored curves on top
plots = get(gca,'Children');
offset = 1;
for i=1:length(r_to_label)
    index = find(rs == r_to_label(i),1);
    if(~isempty(index))
        plots = [plots(end-index+offset); plots];
        plots(end-index+offset) = [];
        offset = offset + 1;
    end
end
set(gca,'Children', plots);

% Add Legend
handle_array = cellfun(@(x) x, handles);
label_array  = cellfun(@(x) ['r = ', n2s(x{1},'%g')], plot_options.labels, 'UniformOutput', false);
if(plot_options.show_legend)
    lg = legend(handle_array,label_array);
    set(lg,'Location','NorthWest','Box','off','Xcolor',[.99 .99 .99],'Ycolor',[.99 .99 .99],'FontWeight','normal');
        % HACK: Push legend up to very top
        p = get(lg,'position'); p(1) = .1 + p(3); p(2) = .8;
        set(lg,'position',p)
end

% Set axis and legend properties
set(gca, 'FontSize',10);
set(gca,'XTick',linspace(-lx,lx,plot_options.num_ticks));
set(gca,'YTick',linspace(-ly,ly,plot_options.num_ticks));
axis(1.05 * plot_options.axis); 

% properly orient figure
view([90 90]); set(gca,'XDir','Reverse');

% set title
if(plot_options.print_title && strcmp(plot_options.type,'accuracy'))
    title(['Accuracy Region $\epsilon = ', n2s(plot_options.epsilon,'%.0e'), '~$ for ',method.type,'SDC$_{', num2str(method.N), '}^{',num2str(method.M),'}$'], 'Interpreter','latex');
elseif(plot_options.print_title && strcmp(plot_options.type,'stability'))
    % extract and format interval angle
    [~,~,rot] = n2s(rs(end),'%g');
    rot = regexprep(rot,'^(-?)1$','$1'); % turn 1 to '' and -1 to '-', dont affect other numbers;    
    rot = [rot, '\cdot'];
    title(['Stability Regions for ',method.type,'SDC$_{', num2str(method.N), '}^{',num2str(method.M),'}$', ' and r $\in ',rot ,'[', num2str(abs(rs(1))), ',' , num2str(abs(rs(end))),']$'], 'Interpreter','latex');
end

% box around plot
set(gca, 'box',plot_options.box);
       
end


function [s,rs,as]=n2s(n,format)
% N2S converts number to string (in a nice looking way).
    rs = ''; as = '';
    if(nargin == 1)
        s = num2str(n,format);
    elseif(format(1) == '%' && (format(end) == 'f' || format(end) == 'g')) % "complex float/compact-float" (remove 0 from purely imaginary numbers, format complex as exponential in latex)
        if(imag(n) == 0)
            s = num2str(n,format);
            rs = num2str(abs(n),format); as = num2str(sign(n));
        elseif(real(n) == 0)
            s = [num2str(imag(n), format),'i'];
            rs = num2str(abs(n),format); as = [num2str(sign(imag(n))), 'i'];
        else
            a = angle(n); r = abs(n);
            [n,m] = rat(a/pi);
            s1 = regexprep(num2str(n),'^(-?)1$','$1'); % turn 1 to '' and -1 to '-', dont affect other numbers; 
            rs = num2str(abs(r),format); as = ['e^{',s1,'\pi i/',num2str(m),'}'];
            s = [rs, as];            
        end
    elseif(format(1) == '%' && format(end) == 'e') % format exponential in latex
        s = [num2str(n,format),'}}'];
        s = strrep(s,'e','\ensuremath{\times 10^{');      
    end
end