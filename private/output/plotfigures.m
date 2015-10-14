function [] = plotfigures(Sample)
% This function creates several useful plot
% keyboard
% Prepare for plotting functions

fprintf(1,'Plotting Figures\n');


%% Plot Tomogram / Image data and movies

% plot histogram of grayscale values
% figure;
histVals = 0.5:1:255.5;
hist(Sample.T(:),histVals);
xlabel( 'Intensity');
ylabel('Counts');
        

% plot 2x2 figure of acf results        
plotfittingresults(Sample);

% plot 3D tomogram with principal planes as surfaces
plottomogramwithplanes(Sample);

% plot principal direction slices
switch Sample.d
    case 3
        plotprincipalplanes(Sample)
    otherwise
end

end


function [hFig] = plotfittingresults(Sample)
        
grayColor = [127,127,127]/256;
scaleFactor = 1/Sample.pxPerMm;

hFig = figure(2);

%%
switch Sample.d
    case 3
        % make 3D slice plot
        xslice = [-(Sample.Ny-1)/2,(Sample.Ny-1)/2]; 
        yslice = [-(Sample.Nx-1)/2,(Sample.Nx-1)/2]; 
        zslice = [-(Sample.Nz-1)/2,(Sample.Nz-1)/2];
        
%         subplot(2,2,1);
        hAx = axes;
        slice(Sample.x1,Sample.x2,Sample.x3,Sample.T,xslice,yslice,zslice,'nearest')
        shading flat;
        axis equal
        axis off
        colormap(Sample.cmapTomogram);
%         title(Sample.name)
        
    case 2
    
%         subplot(2,2,1)
        hAx = axes;
        imagesc(Sample.T);
        hAx(1).YDir = 'normal';
        colormap(gray(256));
        axis equal
        axis off
end

% hTitle(1) = title(sprintf('a) %s',Sample.sampleTitle),'Interpreter','none');
% hTitle(1) = text(0,0,sprintf('a) %s',Sample.sampleTitle),'Interpreter','none');



% plot zeta vs ri
% subplot(2,2,2);
hAx(2) = axes;
plot(hAx(2),Sample.r*scaleFactor,Sample.zeta,'.','Color',grayColor);

% ylim([-0.5,2.5]);
% xlim([0,Sample.maxRadialLag*scaleFactor])
hAx(2).XLim = [0,Sample.maxRadialLag*scaleFactor];

% hTitle(2) = text(0,0,'b) Estimated ACF');
xlabel('r (mm)');
ylabel('\zeta');

% Plot data and best-fit solution in undeformed coordinates
[~,idx] = sort(Sample.r0,'ascend');
hAx(3) = axes;
plot(hAx(3),Sample.r0*scaleFactor,Sample.zeta,'.','Color',grayColor,'DisplayName','observed');
% ylim([-0.5,2.5]);
xlim([0,Sample.maxRadialLag*scaleFactor])

% hTitle(3) = text(0,0,'c) Isotropic ACF');
xlabel('r_0(mm)')
ylabel('\zeta')


% residuals, r0
hAx(4) = axes;        
plot(hAx(4),Sample.r0*scaleFactor,Sample.epsilonHat.unweightedStandardized,'.','Color',grayColor);
        
ylim([-5,5]);
xlim([0,Sample.maxRadialLag*scaleFactor])
hAx(4).YTick = [-5:2.5:5];

% hTitle(4) = text(0,0,'d) Standardized Residuals');
xlabel('r_0(mm)')
ylabel('$\frac{\epsilon}{\sigma_{\epsilon}}$','Interpreter','LaTex')


% arrange axes and titles
PlotOpts = setdefaultplottingopts;
PlotOpts.figureSize = 'fullPage';

hFig = publishfigure(hFig,PlotOpts);
figW = hFig.Position(3); % pos = [left,bot,wid,ht];
figH = hFig.Position(4);

subAxVtSpacing = 3; % cm
subAxHzSpacing = 2; % cm
subAxHeight = 0.5*(figH - 3*subAxVtSpacing);
subAxWidth  = 0.5*(figW - 3*subAxHzSpacing);
subAxHeight = subAxWidth;

hAx(1) = publishfigure(hAx(1),PlotOpts);
hAx(1).Position(1) = subAxHzSpacing;
hAx(1).Position(2) = figH-subAxVtSpacing-subAxHeight;
hAx(1).Position(3) = subAxWidth;
hAx(1).Position(4) = subAxHeight;


hAx(2) = publishfigure(hAx(2),PlotOpts);
hAx(2).Position(1) = subAxWidth + 2*subAxHzSpacing;
hAx(2).Position(2) = figH-subAxVtSpacing-subAxHeight;
hAx(2).Position(3) = subAxWidth;
hAx(2).Position(4) = subAxHeight;

hAx(3) = publishfigure(hAx(3),PlotOpts);
hAx(3).Position(1) = subAxHzSpacing;
hAx(3).Position(2) = figH-2*subAxVtSpacing-2*subAxHeight;
hAx(3).Position(3) = subAxWidth;
hAx(3).Position(4) = subAxHeight;

hAx(4) = publishfigure(hAx(4),PlotOpts);
hAx(4).Position(1) = subAxWidth + 2*subAxHzSpacing;
hAx(4).Position(2) = figH-2*subAxVtSpacing-2*subAxHeight;
hAx(4).Position(3) = subAxWidth;
hAx(4).Position(4) = subAxHeight;

hAx(1).Position(2) = hAx(2).Position(2);



% create new axes to hold annotations
hAxA = axes;
hAxA.Units = 'normalized';
hAxA.Position = [0,0,1,1];
hAxA.Visible = 'off';


titleSpace = 0.5;

hTitle(2) = text(0,0,'b) Estimated ACF');
hTitle(2).Parent = hAxA;
hTitle(2).Units = 'centimeters';
hTitle(2).Position(1) = hAx(2).Position(1);
hTitle(2).Position(2) = hAx(2).Position(2) + subAxHeight + titleSpace;


hTitle(3) = text(0,0,'c) Isotropic ACF');
hTitle(3).Parent = hAxA;
hTitle(3).Units = 'centimeters';
hTitle(3).Position(1) = hAx(3).Position(1);
hTitle(3).Position(2) = hAx(3).Position(2) + subAxHeight + titleSpace;

hTitle(4) = text(0,0,'d) Standardized Residuals');
hTitle(4).Parent = hAxA;
hTitle(4).Units = 'centimeters';
hTitle(4).Position(1) = hAx(4).Position(1);
hTitle(4).Position(2) = hAx(4).Position(2) + subAxHeight + titleSpace;

hTitle(1) = text(0,0,sprintf('a) %s',Sample.sampleTitle));
hTitle(1).Parent = hAxA;
hTitle(1).Units = 'centimeters';
hTitle(1).Position(1) = hTitle(3).Position(1);
hTitle(1).Position(2) = hTitle(2).Position(2);


publishfigure(hTitle(1));
publishfigure(hTitle(2));
publishfigure(hTitle(3));
publishfigure(hTitle(4));



hAx(2).YLabel.Rotation = 0;
hAx(3).YLabel.Rotation = 0;
hAx(4).YLabel.Rotation = 0;

hAx(2).YLabel.FontSize = 24;
hAx(3).YLabel.FontSize = 24;
hAx(4).YLabel.FontSize = 24;

hAx(2).YLabel.Position(1) = -hAx(2).XLim(2)/6;
hAx(3).YLabel.Position(1) = -hAx(3).XLim(2)/6;



end


function [fig] = plottomogramwithplanes(Sample)
switch Sample.d
    case 3
        % Plot 3D Slice with principal directions
        fig = figure;
        hold all

        %... make 3D slice
        xslice = [-(Sample.Ny-1)/2,(Sample.Ny-1)/2]; 
        yslice = [-(Sample.Nx-1)/2,(Sample.Nx-1)/2]; 
        zslice = [-(Sample.Nz-1)/2,(Sample.Nz-1)/2];
        slice(Sample.x1,Sample.x2,Sample.x3,Sample.T,xslice,yslice,zslice,'nearest')
        shading flat;

        % set color map
        colormap(Sample.cmapTomogram);
%         set(gca,'CLim',[0,255])
        

        % plot principal directions
% delete(hLine)        
        hLine(1) = line([0,-2*Sample.vecParSx(1,end)],[0,-2*Sample.vecParSx(2,end)],[0,-2*Sample.vecParSx(3,end)],...
            'LineWidth',2,'Color','black',  'DisplayName','X');
        hLine(2) = line([0,2*Sample.vecParSy(1,end)],[0,2*Sample.vecParSy(2,end)],[0,2*Sample.vecParSy(3,end)],...
            'LineWidth',2,'Color','black','DisplayName','Y');
        hLine(3) = line([0,-2.5*Sample.vecParSz(1,end)],[0,-2.5*Sample.vecParSz(2,end)],[0,-2.5*Sample.vecParSz(3,end)],...
            'LineWidth',2,'Color','black', 'DisplayName','Z');
%         legend(hLine)
        
% delete(hT);
        hT(1) = text(-2.3*Sample.vecParSx(1,end),-2.1*Sample.vecParSx(2,end),-2.2*Sample.vecParSx(3,end),'X','FontSize',18);
        hT(2) = text(2.1*Sample.vecParSy(1,end),2*Sample.vecParSy(2,end),2.3*Sample.vecParSy(3,end),'Y','FontSize',18);
        hT(3) = text(-2.5*Sample.vecParSz(1,end),-2.5*Sample.vecParSz(2,end),-2*Sample.vecParSz(3,end),'Z','FontSize',18);

        % plot principal planes
        hSurf = surf(Sample.SzNorm.x.deformed,Sample.SzNorm.y.deformed,Sample.SzNorm.z.deformed);
        set(hSurf,'EdgeAlpha',0.3,'FaceAlpha',0.3,'Alphadata',0.3 )

        hSurf =surf(Sample.SxNorm.x.deformed,Sample.SxNorm.y.deformed,Sample.SxNorm.z.deformed);
        set(hSurf,'EdgeAlpha',0.3,'FaceAlpha',0.3,'Alphadata',0.3 )

        hSurf =surf(Sample.SyNorm.x.deformed,Sample.SyNorm.y.deformed,Sample.SyNorm.z.deformed);
        set(hSurf,'EdgeAlpha',0.3,'FaceAlpha',0.3,'Alphadata',0.3 )

        % figure props
        axis equal
%         xlabel('x')
%         ylabel('y')
%         zlabel('z')

        axis off;
        set(gca,'View',[-42,30]) % -57,40 for 90727-3

    case 2
        %... plot 2D image
        fig = figure; 
        pcolor(Sample.T); shading flat;
        axis equal
        axis off
        colormap(gray(256))
        scaleBarLengthPixels = 100;    
        line([1,scaleBarLengthPixels],[5,5],[0,0],'Color','Black')
end

end


function [] = plotprincipalplanes(Sample)
% function to plot each principal slice of tomogram in both deformed and
% undeformed coordinates

PlotOpts = setdefaultplottingopts;

%%
% Tomogram, deformed coordinates
hFig(1) = figure;
hAx(1) = axes;
pcolor(Sample.SyNorm.plane.tomogram.deformed);
shading flat;
axis equal; 
axis off;
title(sprintf('X-Z Section\n Deformed Coordinates'));
colormap(Sample.cmapTomogram)
set(gca,'CLim',[0,255])


hFig(2) = figure;
hAx(2) = axes;

pcolor(Sample.SzNorm.plane.tomogram.deformed);
shading flat;
axis equal;
axis off;
title(sprintf('X-Y Section\nDeformed Coordinates'));
colormap(Sample.cmapTomogram)
set(gca,'CLim',[0,255])



hFig(3) = figure;
hAx(3) = axes;
pcolor(Sample.SxNorm.plane.tomogram.deformed);
shading flat;
axis equal;
axis off;
title(sprintf('Y-Z Section\nDeformed Coordinates'));
colormap(Sample.cmapTomogram)
set(gca,'CLim',[0,255])

% Tomogram, undeformed coordinates
hFig(4) = figure;
hAx(4) = axes;

pcolor(Sample.SyNorm.plane.tomogram.undeformed);
shading flat;
axis equal; 
axis off;
title(sprintf('X-Z Section\nUndeformed Coordinates'));
colormap(Sample.cmapTomogram)
set(gca,'CLim',[0,255])

hFig(5) = figure;
hAx(5) = axes;
pcolor(Sample.SzNorm.plane.tomogram.undeformed);
shading flat;
axis equal;
axis off;
title(sprintf('X-Y Section\nUndeformed Coordinates'));
colormap(Sample.cmapTomogram)
set(gca,'CLim',[0,255])

hFig(6) = figure;
hAx(6) = axes;
pcolor(Sample.SxNorm.plane.tomogram.undeformed);
shading flat;
axis equal;
axis off;
title(sprintf('Y-Z Section\nUndeformed Coordinates'));
colormap(Sample.cmapTomogram)
set(gca,'CLim',[0,255])

%% loop through figures and add labels
for iF = 1:6
    
    ax(iF) = axes;
    ax(iF).Parent = hFig(iF);
    publishfigure(hFig(iF));
    publishfigure(hAx(iF));

    ax(iF).Units = 'centimeters';
    plot(ax(iF),[0,0],[0,1],'-k');
    hold(ax(iF),'on');
    plot(ax(iF),[0,1],[0,0],'-k');

    pos = hAx(1).Position;
    % axis(ax,'equal');
    ax(iF).Position(1) = pos(1);
    ax(iF).Position(2) = pos(2);
    ax(iF).Position(3) = 0.25*pos(3);
    ax(iF).Position(4) = 0.25*pos(3);
    ax(iF).Visible = 'off';

    txtZ(iF) = text(0,1.2,'Z');
    txtZ(iF).FontSize = 14;
    txtZ(iF).HorizontalAlignment = 'center';
    txtZ(iF).Parent = ax(iF);

    txtX(iF) = text(1.2,0,'X');
    txtX(iF).FontSize = 14;
    txtX(iF).HorizontalAlignment = 'center';
    txtX(iF).Parent = ax(iF);

end

txtZ(2).String = 'Y';

txtZ(3).String = 'Z';
txtX(3).String = 'Y';

txtX(5).String = 'X';
txtZ(5).String = 'Y';

txtX(6).String = 'Y';
txtZ(6).String = 'Z';

%% ACF (rho)
% levels = 5;
% cntrColor = 'Black';
% figure;
% subplot(3,2,1);
% pcolor(Sample.SyNorm.plane.rho.deformed);
% hold all;
% contour(Sample.SyNorm.plane.rho.deformed,levels,'LineColor',cntrColor)
% shading flat;
% axis equal; 
% axis off;
% title('X-Z Section, Deformed Coordinates')
% colormap(Sample.cmapRho)
% zoom(4)
% 
% subplot(3,2,3);
% pcolor(Sample.SzNorm.plane.rho.deformed);
% hold all;
% contour(Sample.SzNorm.plane.rho.deformed,levels,'LineColor',cntrColor)
% shading flat;
% axis equal;
% axis off;
% title('X-Y Section, Deformed Coordinates')
% colormap(Sample.cmapRho)
% zoom(4)
% 
% subplot(3,2,5);
% pcolor(Sample.SxNorm.plane.rho.deformed);
% hold all;
% contour(Sample.SxNorm.plane.rho.deformed,levels,'LineColor',cntrColor)
% shading flat;
% axis equal;
% axis off;
% title('Y-Z Section, Deformed Coordinates')
% colormap(Sample.cmapRho)
% zoom(4)
% 
% % Rho, undeformed coordinates
% % figure;
% subplot(3,2,2);
% pcolor(Sample.SyNorm.plane.rho.undeformed);
% hold all;
% contour(Sample.SyNorm.plane.rho.undeformed,levels,'LineColor',cntrColor);
% shading flat;
% axis equal; 
% axis off;
% title('X-Z Section, Undeformed Coordinates')
% colormap(Sample.cmapRho)
% zoom(4)
% 
% subplot(3,2,4);
% pcolor(Sample.SzNorm.plane.rho.undeformed);
% hold all;
% contour(Sample.SzNorm.plane.rho.undeformed,levels,'LineColor',cntrColor);
% shading flat;
% axis equal;
% axis off;
% title('X-Y Section, Undeformed Coordinates')
% colormap(Sample.cmapRho)
% zoom(4)
% 
% subplot(3,2,6);
% pcolor(Sample.SxNorm.plane.rho.undeformed);
% hold all;
% contour(Sample.SxNorm.plane.rho.undeformed,levels,'LineColor',cntrColor);
% shading flat;
% axis equal;
% axis off;
% title('Y-Z Section, Deformed Coordinates')
% colormap(Sample.cmapRho)
% zoom(4)


end



