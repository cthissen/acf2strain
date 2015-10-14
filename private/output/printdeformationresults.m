function [] = printdeformationresults(Sample)

defParams = Sample.defParams;
fid = Sample.fid;

% fid = 1;
dualfprintf(fid,sprintf('\n======================GOODNESS OF FIT TESTS=======================\n'));
dualfprintf(fid,sprintf(' --Fraction of total variance accounted for by model fit--\n'));
dualfprintf(fid,sprintf('R^2 Value: %3.2f\n',Sample.R2));
dualfprintf(fid,sprintf(' --Durbin-Watson Test for serial correlation--\n'));
dualfprintf(fid,sprintf('D-W Statistic: %3.2f\n',Sample.dw));
% dualfprintf(fid,sprintf('D-W Statistic: %3.2f    p-value for null hypothesis: %3.2f \n',Sample.dw,NaN));
    


dualfprintf(fid,sprintf('\n\n===========BEST-FIT DEVIATORIC STRAIN ELLIPSOID (normalized to Sv=1)==========\n'));
dualfprintf(fid,sprintf('Axis \t -Stretch- \t ------- Stretch Standard Errors -------    LnStrain',Sample.pValue*100));
dualfprintf(fid,sprintf('\n  \t  (Lf/Li) \t  Lw Lmt  Hi Lmt    -SE    +SE   Ave SE       (cNp) \n'));
%... print results for each axis:
for i = 1:3
    vector = defParams.vector{i};
dualfprintf(fid,sprintf('%i  \t  %3.2f \t\t  %3.2f     %3.2f     %3.2f   %3.2f   %3.2f       %+3.2f     \n',...
    i,...
    vector.stretch.mag, vector.stretch.min,vector.stretch.max,...
    vector.stretch.negSE, vector.stretch.posSE, vector.stretch.aveSE,...
    vector.mag*100));
end

% print isotropy test results
dualfprintf(fid,sprintf('\n\n===============================ISOTROPY TESTS================================\n'));
dualfprintf(fid,sprintf(' --Probability that difference in principal stretches is due to chance alone--\n'));
dualfprintf(fid,sprintf('\t\tX-Z: %3.0f \t X-Y: %3.0f \t Y-Z: %3.0f\n',...
    defParams.isoTest(1,3), defParams.isoTest(1,2), defParams.isoTest(2,3)));
   
    
% print directional results
dualfprintf(fid,sprintf('\n\n============================PRINCIPAL DIRECTIONS=============================\n'));
dualfprintf(fid,sprintf(' --Principal Axis--     ------Error Ellipse in Degrees (%2.0f percent CI)------\n',Sample.pValue*100));
dualfprintf(fid,sprintf('   Axis   Trd, Pl       Plane   Trd, Pl  Radius    Trd, Pl  Radius    Twist\n'));
%... print results for each axis:
for i = 1:3
    vector = defParams.vector{i};
    otherAxis = setdiff(1:3,i);
    

dualfprintf(fid,sprintf('   %i:    %+04.0f, %+3.0f        %i%i:  %+04.0f, %+3.0f  %04.1f    %+04.0f, %+3.0f  %04.1f     %+3.2f\n',...
    i,...
    vector.trendDegrees,              vector.plungeDegrees,...
    otherAxis(1),otherAxis(2),...
    vector.majorEllipse.trendDegrees, vector.majorEllipse.plungeDegrees, vector.majorEllipse.magDegrees,...
    vector.minorEllipse.trendDegrees, vector.minorEllipse.plungeDegrees, vector.minorEllipse.magDegrees,...
    vector.twistDegrees));
end

%% Create Latex File
% natural strains
Ex = Sample.defParams.vector{1}.mag;
Ey = Sample.defParams.vector{2}.mag;
Ez = Sample.defParams.vector{3}.mag;
Ev = Ex + Ey + Ez;

errEx = Sample.defParams.vector{1}.ci;
errEy = Sample.defParams.vector{2}.ci;
errEz = Sample.defParams.vector{3}.ci;
errEv = sqrt(errEx^2 + errEy^2 + errEz^2);

% deviatoric stretches
Sx = exp(Ex);
Sy = exp(Ey);
Sz = exp(Ez);
Sv = exp(Ev);

% calculate average CI's for stretches
SxMax =  exp(Ex + errEx);
SyMax =  exp(Ey + errEy);
SzMax =  exp(Ez + errEz);
SvMax =  exp(Ev + errEv);

SxMin =  exp(Ex - errEx);
SyMin =  exp(Ey - errEy);
SzMin =  exp(Ez - errEz);
SvMin =  exp(Ev - errEv);

SxnegCI = Sx - SxMin;
SynegCI = Sy - SyMin;
SznegCI = Sz - SzMin;
SvnegCI = Sv - SvMin;


SxposCI = SxMax - Sx;
SyposCI = SyMax - Sy;
SzposCI = SzMax - Sz;
SvposCI = SvMax - Sv;


SxaveCI = sqrt((SxposCI^2 + SxnegCI^2)/2);
SyaveCI = sqrt((SyposCI^2 + SynegCI^2)/2);
SzaveCI = sqrt((SzposCI^2 + SznegCI^2)/2);
SvaveCI = sqrt((SvposCI^2 + SvnegCI^2)/2);

fid=Sample.fid;
i=1;
fprintf(fid,'\n\nLatex-formatted Table:\n');
fprintf(fid,['& \\textit{X} & $%03.2f (%2.0f)$ & $%+06.2f (%2.0f)$ &$%03.0f/%02.0f$',...
             '&$%1.1f$ &$%1.1f$ &$%02.0f$ &Y-Z: $%1.0f$&\\multirow{4}{*}{%03.0f}',...
             '&\\multirow{4}{*}{%04.2f} &\\multirow{4}{*}{%04.2f} \\\\\n'],...
             Sx, SxaveCI, Ex*100, errEx*10000,...
			 Sample.defParams.vector{i}.trendDegrees,...
			 Sample.defParams.vector{i}.plungeDegrees,...
             Sample.defParams.vector{i}.majorEllipse.magDegrees,...
             Sample.defParams.vector{i}.minorEllipse.magDegrees,...   
             Sample.defParams.vector{i}.twistDegrees,...
             Sample.defParams.isoTest(2,3), Sample.Rs_micron,Sample.R2,...
             Sample.dw);
i=2;         
fprintf(fid,['& \\textit{Y} & $%03.2f (%2.0f)$ & $%+06.2f (%2.0f)$ &$%03.0f/%02.0f$',...
             '&$%1.1f$ &$%1.1f$ &$%02.0f$ &X-Z: $%1.0f$ \\\\\n'],...
             Sy, SyaveCI, Ey*100, errEy*10000,...         
			 Sample.defParams.vector{i}.trendDegrees,...
			 Sample.defParams.vector{i}.plungeDegrees,...
			 Sample.defParams.vector{i}.majorEllipse.magDegrees,...
             Sample.defParams.vector{i}.minorEllipse.magDegrees,...             
             Sample.defParams.vector{i}.twistDegrees,...
             Sample.defParams.isoTest(1,3));
i=3;        
fprintf(fid,['& \\textit{Z} & $%03.2f (%2.0f)$ & $%+06.2f (%2.0f)$ &$%03.0f/%02.0f$',...
             '&$%1.1f$ &$%1.1f$ &$%02.0f$ &X-Y: $%1.0f$ \\\\\n'],...
             Sz, SzaveCI, Ez*100, errEz*10000,... 
			 Sample.defParams.vector{i}.trendDegrees,...
			 Sample.defParams.vector{i}.plungeDegrees,...                          
             Sample.defParams.vector{i}.majorEllipse.magDegrees,...
             Sample.defParams.vector{i}.minorEllipse.magDegrees,...             
             Sample.defParams.vector{i}.twistDegrees,...
             Sample.defParams.isoTest(1,2));         

fprintf(fid,'& \\textit{V} & $%03.2f(%2.0f)$ & $%+06.2f(%2.0f)$ &--- &--- &--- &--- &--- \\\\\n',...
            Sv, SvaveCI, Ev*100, errEv*10000);

end