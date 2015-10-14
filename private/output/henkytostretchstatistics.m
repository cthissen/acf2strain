function [stat] = henkytostretchstatistics(stat)
% parse output variable from hextstatistics function into confidence
% intervals appropriate for estimation of the Hencky tensor
% adds stretch structure to stat.vector


% parse Hext analysis results
for i = 1:3
    % standard errors for stretch
    stat.vector{i}.stretch.mag = exp(stat.vector{i}.mag);
    
    stat.vector{i}.stretch.max = exp(stat.vector{i}.mag + stat.vector{i}.se);
    stat.vector{i}.stretch.min = exp(stat.vector{i}.mag - stat.vector{i}.se);
    
    stat.vector{i}.stretch.negSE = stat.vector{i}.stretch.mag - stat.vector{i}.stretch.min;
    stat.vector{i}.stretch.posSE = stat.vector{i}.stretch.max - stat.vector{i}.stretch.mag;
    stat.vector{i}.stretch.aveSE = sqrt((stat.vector{i}.stretch.posSE^2 + stat.vector{i}.stretch.negSE^2)/2);

    

end


end