function [] = savefigure_cjt(hFig,figName,varargin)
%[] = saveFigure(hFig,figName,varargin)
% save figures using export_fig

%% check inputs
narginchk(2,Inf);
validateattributes(figName,{'char'},{});

if nargin == 2
    % default is save current figure as png, with transparent background
    imgType = '-png';

else
    imgType = varargin{1};
end

%%
switch imgType
    case '-png'
        %... get al subplots
%         ax = findall(hFig,'Type','Axes');
%         for i = 1:numel(ax)            
%             set(ax(i),'Color','none');
%         end            
        set(hFig,'Color','white');
        export_fig(hFig,figName,'-png','-transparent','-m3');

    case '-pdf'
        
        set(hFig,'Color','white');
        export_fig(hFig,figName,'-pdf');
end
    







end