function [PlotOpts] = setdefaultplottingopts()
% [PlotOpts] = setdefaultplottingopts()

% figure options
PlotOpts.f                  = 1; % figure number
PlotOpts.hFig               = []; % figure handle
PlotOpts.totalWidth_cm      = 19;   % total figure width in cm (portrait full page figure = 19x23cm)
PlotOpts.totalHeight_cm     = 24;
PlotOpts.figureSize         = 'oneColumn';

% axes options
% PlotOpts.ax                 = axes; % new axes
PlotOpts.fontSize           = 14;
PlotOpts.fontName           = 'Helvetica';
PlotOpts.axLineWidth        = 1;

% line options
PlotOpts.lineWidth          = 1.5;

% scatterplot options
PlotOpts.scatterPointSize   = 50;    % default scaling for point size in q-plot


% PlotOpts.subaxesSpacing_cm  = 0.25; % spacing between stereograms (cm)
% PlotOpts.colorBarSpacing_cm = 0.75; % spacing at bottom of figure for colorbar (cm)
% PlotOpts.probabilityCalc    = 'parameterized'; % parameterized or bootstrap probability


end