function [] = makemovie(Sample)
    %% Movie Capture
    
    Folder = loadDirectoryInfo;
    l = Sample.Nz;
    writerObj = VideoWriter([Folder.results,Sample.name]);
    open(writerObj);
    hFig = figure(30);
    axis equal
    set(hFig,'Renderer','zbuffer');
    set(gca,'nextplot','replacechildren');
    set(gca,'YDir','reverse');
    
    cLims(1) = min(255*Sample.T(:));
    cLims(2) = max(255*Sample.T(:));
    
    for i = 1:l
        imagesc(255*squeeze(Sample.T(:,:,i)))
            axis equal
            axis off;
            colormap(Sample.cmapTomogram)
            set(gca,'CLim',cLims)
        frame = getframe;
        writeVideo(writerObj,frame);
    end
    close(writerObj);
    
    close(hFig);
end



