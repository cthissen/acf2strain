function [h] = publishfigure(h,varargin)
% [hFig] = publishfigure(hFig,varargin)
% this function sets some basic figure properties appropriate for
% submitting for publications. 

%% check inputs
narginchk(1,2); 

if nargin == 1
    PlotOpts = setdefaultplottingopts;
else
    PlotOpts = varargin{1};
end
handleType = h(1).Type;
switch handleType
    case 'figure'
        % first set size
        figureSize = PlotOpts.figureSize;
        h.Units = 'centimeters';
        h.Color = 'w';

        switch figureSize
            case 'oneColumn'
                position = [0,0,9,9];        
            case 'halfPage'
                position = [0,0,14,24];        
            case 'fullPage'
                position = [0,0,19,24];
        end
        h.Position = position;
        movegui(h,'northwest');
    case 'axes'

        % set font name and font size
        %... select all axes in figure
        h.FontName  = PlotOpts.fontName;
        h.FontSize  = PlotOpts.fontSize;
        h.LineWidth = PlotOpts.axLineWidth;
        h.Units = 'centimeters';
        
    case 'colorbar'
        h.FontName = PlotOpts.fontName;
        h.FontSize = PlotOpts.fontSize;
        
    case 'line'
        for ih = 1:numel(h)
            h(ih).LineWidth = PlotOpts.lineWidth;
        end
        
    case 'text'
        h.FontName  = PlotOpts.fontName;
        h.FontSize  = PlotOpts.fontSize;  
    otherwise
        warning('publishFigure: No options set for this type of handle object');
end

end