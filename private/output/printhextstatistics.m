function [] = printhextstatistics(stat)
fid = 1;
%%

% prints Hext results to screen, as output by parsehextstatistics
%... modeled after Hext, table 12.3
fprintf(fid,sprintf('\n\n=============== Principal Axes, Hext Table 12.3 ===============\n'));
fprintf(fid,sprintf('Axis \t -Length-    Direction cosines     sigma    %02.0f per. CI  \n',stat.pvalue*100));
%... print results for each axis:
for i = 1:3
    vector = stat.vector{i};
fprintf(fid,sprintf('%i  \t  %3.2f   %+3.3f %+3.3f %+3.3f    %3.2f   %3.2f %3.2f    \n',...
    i,...
    vector.mag, ...
    vector.dirCosines,...
    vector.se, vector.mag-vector.ci, vector.mag+vector.ci));
end


% print isotropy test results
fprintf(fid,sprintf('\n\n====================ISOTROPY TESTS Hext Table 12.2?===================\n'));
fprintf(fid,sprintf('-Probability that difference in principal axes is due to chance alone-\n'));
fprintf(fid,sprintf('\t\t1-3: %3.2f \t 1-2: %3.2f \t 2-3: %3.2f\n',...
    stat.isoTest(1,3), stat.isoTest(1,2), stat.isoTest(2,3)));
   
%%    
% print directional results
fprintf(fid,sprintf('\n\n============================PRINCIPAL DIRECTIONS=============================\n'));
fprintf(fid,sprintf('--Principal Axis--  --------------Error Ellipse in Degrees (%3.0f per. CI)-------\n',stat.pvalue*100));
fprintf(fid,sprintf('   Axis   Trd, Pl       Plane   Trd, Pl  Radius    Trd, Pl  Radius    Twist\n'));
%... print results for each axis:
for i = 1:3
    vector = stat.vector{i};
    otherAxis = setdiff(1:3,i);

fprintf(fid,sprintf('   %i:    %+04.0f, %+3.0f        %i%i:  %+04.0f, %+3.0f  %04.1f    %+04.0f, %+3.0f  %04.1f     %+3.2f\n',...
    i,...
    vector.trendDegrees,              vector.plungeDegrees,...
    otherAxis(1),otherAxis(2),...
    vector.majorEllipse.trendDegrees, vector.majorEllipse.plungeDegrees, vector.majorEllipse.magDegrees,...
    vector.minorEllipse.trendDegrees, vector.minorEllipse.plungeDegrees, vector.minorEllipse.magDegrees,...
    vector.twistDegrees));
end



end