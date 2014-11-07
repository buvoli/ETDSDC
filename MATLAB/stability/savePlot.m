function savePlot(file_name, save_options)
%SAVEPLOT Prints file as specified by save_options
% Useful Documentation
% http://www.mathworks.com/help/matlab/ref/figure_props.html
% http://www.mathworks.com/help/matlab/ref/print.html

set(gcf,'renderer',save_options.renderer);
set(gcf,'PaperUnits',save_options.units);
set(gcf, 'PaperPosition', save_options.position);
set(gcf, 'PaperPositionMode', 'manual'); % force size to PaperPosition specification

print('-depsc2', sprintf('-r%d', save_options.dpi), strcat(file_name, '.eps'));
end