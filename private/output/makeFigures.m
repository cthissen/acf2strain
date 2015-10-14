function [] = makeFigures(Sample,Folder)


    if max(Sample.T(:)) == 1
        Sample.T = 255*Sample.T;
    end
    switch Sample.flag.plotFigures
        case 'false'
            % make no plots

        case 'true'
            % prepare for plotting functions
            %... get planes with axes parallel to principal direction
            Sample = getprincipalplanes(Sample);

            %... get cmap for whole tomogram volume
            cmap = colormap(gray(256));
            Sample.cmapTomogram = cmapscale((Sample.T),cmap,0.5);

            %... get cmap for sensitivity plot and rho plots
            cmap = coolwarm(256);
            Sample.cmapRho = cmapscale(Sample.rho,      cmap,0.5);

            % plot figures
            plotfigures(Sample)

            % save figures
            switch Sample.flag.saveFigures
                case 'true'
                    PlotOpts = setdefaultplottingopts;
                    h = get(0,'children');                    
                    switch Sample.dimensions

                        case 2
                            publishfigure(h(1),PlotOpts);
                            publishfigure(h(3),PlotOpts);
                            PlotOpts.figureSize = 'fullPage';
                            publishfigure(h(2),PlotOpts);
                            
                            figFileName = sprintf('%s%s_%02i',Folder.results,Sample.name,1);
                            savefigure_cjt(h(1),figFileName,'-png');
                            
                            figFileName = sprintf('%s%s_%02i',Folder.results,Sample.name,2);
%                             savefigure_cjt(h(2),figFileName,'-png');
                            set(h(2),'Color','white');
                            export_fig(h(2),figFileName,'-png');
                            
                            figFileName = sprintf('%s%s_%02i',Folder.results,Sample.name,3);
                            savefigure_cjt(h(3),figFileName,'-png');                            
                            
                        case 3
                            h = get(0,'children');
                            for i = 1:length(h)
                                figFileName = sprintf('%s%s_%02i',Folder.results,Sample.name,i);
                                savefigure_cjt(h(i),figFileName,'-png');
                            end

                            % workaround: do fig 8 separate
                            figFileName = sprintf('%s%s_%02i',Folder.results,Sample.name,8);
                            export_fig(h(8),figFileName,'-png');                            
                    end
                                
                   
                    
                case 'false'

                otherwise
                    error('Specify true or false for Sample.flag.saveFigures');
            end

    otherwise
        error('Specify true or false for Flag.plotFigures');
    end

end