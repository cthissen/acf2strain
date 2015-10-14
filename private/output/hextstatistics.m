function [stat] = hextstatistics(Emat,Cov,df,pvalue)
% Using estimates for the components of a 3x3 tensor, find the confidence
% intervals for the principal components (e.g. stretches) and confidence
% intervals for the directions. Also includes an isotropy test.

% Need to specify 6x6 covariance matrix, the degrees of freedom, and a
% pvalue for the confidence intervals. Hext method assumes tensor
% components are normally distributed.


% V0.91 CJT, 11/6/14
    % changed the calculation from using the L matrix to using the direct
    % equations in Hext. Also added loops so output is given as vector{i}
    % structure for i=1:3
    % added aoperator subfunction

% V0.9 CJT, 11/4/14
    % Matlab implementation of Hext Statistics. Originally programmed as
    % part of STRAIN3D.BAS. 
    % 7/1/14 - moved to StrainACF folder. 
        % modified to accept df as input parameter
        % note that Binv does not need to be the Finger Tensor (inverse
        % Cauchy-Green tensor)
    % 10/20/14 modified to accept p-value (probability level for confidence
    % intervals)
    % 11/4/14 modified to use more straightforward variable names
    
% some use references:
    % Hext, 1963
    % Jelinek, 1977, 1978
        % 1978 might be appropriate for "regional average strain"
    % Constable and Tauxe, 1990
    % Lienert, 1991
    % Owens, 2000 GJI 142, 527-538
        
% Example Input Parameters    
% Binv = diag([2^2,0.5^2,1^2]);
% Cov = diag([0.1,0.1,0.1,0.1,0.1,0.1]);
% N = 4; % in strain3D this is the number of section ellipses. 



% Set constants
RAD = pi/180;
TINY = 1e-15; % double precision limit


%... Find principal values and directions for Binv
[T,E] = sortEig(Emat,'descend'); % 

%... Calculate the trend and plunge of the eigenvectors
for j = 1:3
    % ensure that all eigenvectors are lower hemisphere
    if T(3,j) < 0
        T(1,j) = -T(1,j);
        T(2,j) = -T(2,j);
        T(3,j) = -T(3,j);
    end
    [trend(j),plunge(j)] = xyz2tp([T(1,j),T(2,j),T(3,j)]);
    
    vector{j}.mag = E(j,j);
    vector{j}.trendDegrees  = adjang(trend(j),0,2*pi)/RAD;
    vector{j}.plungeDegrees = plunge(j)/RAD;

end
trend  = adjang(trend, 0,2*pi);

%% Calculate Errors for Principal Stretches (values)
%... Calculate the +/- SE for principal stretch, corrected for small sample

TwoQ = 1 - pvalue;
% finds the root ofTwoQ - B(0.5V,0.5, (1V/(1V+0.1^2))) where B is the
% incomplete beta function (the CDF for a binomial is the incomplete beta
% function)
DF = df;
factor = Tpoints(TwoQ,DF);
for i = 1:3 % loop through eigenvalues
    iGamma = T(:,i);
    aii = aoperator(iGamma, iGamma);
    
    EstandardErrors = factor * sqrt(aii'*Cov*aii); % eq 5.5. 
    vector{i}.se = sqrt(aii'*Cov*aii);
    vector{i}.ci = EstandardErrors;
end

%% Calculate Errors for Principal Directions
% Calculate the error ellipses for principal directions.
% Corrected for small-sample bias using F ratio. 
factor = Fpoints(pvalue,2,DF);
for i = 1:3 % loop through eigenvectors. Recall they descend
    % get indices
    %... j is always given largest principal direction that is not iGamma    
    if i == 1
        j = 2;
        k = 3;
    elseif i == 2
        j = 1;
        k = 3;
    elseif i == 3
        j = 1; 
        k = 2;
    end
    

        
    % construct Ediff
    Ediffij = 1/(E(i,i) - E(j,j));
    Ediffik = 1/(E(i,i) - E(k,k));
    
    % construct gamma vectors
    iGamma = T(:,i);
    jGamma = T(:,j);
    kGamma = T(:,k);
    
    % construct a vectors
    aij = aoperator(iGamma,jGamma);
    aik = aoperator(iGamma,kGamma);
    
    % Eq 4.17
    RHS = [Ediffij*aij ,  Ediffik*aik];
    
    LHS = [Ediffij*aij' ;...
           Ediffik*aik'];
           
    Wi = LHS * Cov * RHS;
    
    
    %... check for roundoff error
    Wi(Wi<TINY) = 0;
    
    % confidence intervals are related to the eigenvectors and eigenvalues of Wi
    %... sorted as j,k. (Previously, was sorted according to magnitude of
    % error, at expense of loss of information about which axis)
    [c,lami] = sortEig(Wi,'descend');
    
    % angular distance of ellipse axes
    %... why factor of 2?
    jDirR = atan2(sqrt(2*factor*lami(1,1)),1)/RAD;
    kDirR = atan2(sqrt(2*factor*lami(2,2)),1)/RAD;
    
    % calculate directions of error ellipses in plane perpendicular to iGamma
    %... twist is rotation of major axis from j direction in plane.
    twist = atan(c(2,1)/c(1,1)); % point inversion symmetry, 180 deg range is fine
    
    %... to compute direction of ellipse axis in 3D coordinates, rotate the j or k
    % direction around the i direction by twist degrees
    %...rotate = rotate(rotationaxis, degrees, vector to rotate), radians
    [jDirTrend,jDirPlunge] = rotate_yale(trend(i),plunge(i),twist,trend(j),plunge(j));
    [kDirTrend,kDirPlunge] = rotate_yale(trend(i),plunge(i),twist,trend(k),plunge(k));
    
    % save for output
    vector{i}.majorEllipse.magDegrees    = jDirR;
    vector{i}.majorEllipse.trendDegrees  = adjang(jDirTrend,0,2*pi)/RAD;
    vector{i}.majorEllipse.plungeDegrees = jDirPlunge/RAD;

    vector{i}.minorEllipse.magDegrees    = kDirR;
    vector{i}.minorEllipse.trendDegrees  = adjang(kDirTrend,0,2*pi)/RAD;
    vector{i}.minorEllipse.plungeDegrees = kDirPlunge/RAD;

    vector{i}.twistDegrees = twist/RAD;
   
    
    
end

%% Isotropy test
% Hext eq 6.2
gamma1 = T(:,1);
gamma2 = T(:,2);
gamma3 = T(:,3);

a11 = aoperator(gamma1,gamma1);
a22 = aoperator(gamma2,gamma2);
a33 = aoperator(gamma3,gamma3);

diffStat12 =(E(1,1)-E(2,2)) / sqrt((a11-a22)' * Cov * (a11-a22));
TwoQ = betainc(df/(df+diffStat12^2),0.5*df,0.5);
isoTest(1,2) = 100*TwoQ;

diffStat13 =(E(1,1)-E(3,3)) / sqrt((a11-a33)' * Cov * (a11-a33));
TwoQ = betainc(df/(df+diffStat13^2),0.5*df,0.5);
isoTest(1,3) = 100*TwoQ;

diffStat23 =(E(2,2)-E(3,3)) / sqrt((a22-a33)' * Cov * (a22-a33));
TwoQ = betainc(df/(df+diffStat23^2),0.5*df,0.5);
isoTest(2,3) = 100*TwoQ;
 
%% Parse Statistical Results Print out results
vector{1}.dirCosines = T(:,1);
vector{2}.dirCosines = T(:,2);
vector{3}.dirCosines = T(:,3);


stat.vector = vector;
stat.isoTest = isoTest;
stat.pvalue = pvalue;
stat.df = df;

end

%% subfunctions 
function [aij] = aoperator(jVec, hVec)
% construct 6x1 column vector a. See definiont under Hext eq 3.6

    aij = [jVec(1)*hVec(1);...
           jVec(2)*hVec(2);...
           jVec(3)*hVec(3);...
           jVec(2)*hVec(3) + jVec(3)*hVec(2);...
           jVec(3)*hVec(1) + jVec(1)*hVec(3);...
           jVec(1)*hVec(2) + jVec(2)*hVec(1)];
           

end

%% BASIC Functions from Yale Deformation that have been replaced with MATLAB Versions

function factor = Tpoints(TwoQ,V)
% '=============================================================================
% '... Using bisection, find the T value, called Tpoints, at a specified
% ' two-tail signficance, called TwoQ, with V% degrees of freedom.
% ' The result is refined until the accuracy < Yacc. The bracket used here,
% ' {1,100} is good for values of TwoQ from 100 to 1%.
% ' Modified from RTBIS in Press et al., 1986, p. 247.
% '=============================================================================
% Comparing the matlab equation for the incomplete beta function with the
% implementation used in the STAT.SUB, 
% it is clear that A = Z, B = W, X = X
% Matlab:  betainc(X,Z,W). Be sure to reconvert when using xmin
% Stat.sub: Betai (A, B, X)
    % BT = EXP(GammLn(A + B) - GammLn(A) - GammLn(B) + A * LOG(X) + B * LOG(1! - X))
    % Betai(.5 * V%, .5, 1! * V% / (1! * V% + Y1 * Y1))



fnc = @(x) TwoQ - betainc(x,0.5*V,0.5);
y1 = 0.1;
y2 = 100;


x01 = 1*V/(1*V+y1^2);
x02 = 1*V/(1*V+y2^2);
x = fzero(fnc,[x01,x02]);
% convert x to factor 
factor = sqrt(V/x-V);

end

function [factor] = Fpoints(Prob,v1,v2)
% '=============================================================================
% '... Using bisection, find the F value, called Fpoints, at a specified
% ' cumulative probability, called Prob. F is equal to a ratio of variances
% ' with the numerator and denominator having V1% and V2% degrees of freedom,
% ' respectively. The result is refined until the accuracy < Yacc. The
% ' bracket used here {1,100} is good for probabilities between 1% and 99%,
% ' as long the V1% and V2% remain greater than about 1 or 2.
% ' Modified from RTBIS in Press et al., 1986, p. 247.
% '=============================================================================

% Comparing the matlab equation for the incomplete beta function with the
% implementation used in STAT.SUB, 
% it is clear that A = Z, B = W, X = X
% Matlab:  betainc(X,Z,W) (betainc(X,A,B)
% Stat.sub: Betai (A, B, X)
    % BT = EXP(GammLn(A + B) - GammLn(A) - GammLn(B) + A * LOG(X) + B * LOG(1! - X))

%%
V1 = v1;
V2 = v2;
    
y1 = 1;
y2 = 50;
p = Prob;

if p < 0.5
    p = 1-p;
    tmp = v1;
    v1 = v2;
    v2 = tmp;
end


%... nested function to minimize
    function y = myfnc(x,v)
        v1 = v(1);
        v2 = v(2);
        tmp = betainc(x,0.5*v2,0.5*v1);
        if tmp > 0.5
            tmp = 1-tmp;
        end
        y = 1-p-tmp;
    end
    v = [v1 v2];
    fnc = @(x)myfnc(x,v); % now make function only a function of x
           
%... attempt to bracket zero
x01 = v2/(v2+v1*y1);
x02 = v2/(v2+v2*y2);

y = fzero(fnc,[x01,x02]);
% convert y 
y = (1/v1)*(v2/y-v2);

if Prob < 0.5
    y = 1/y;
    tmp = v1;
    v1 = v2;
    v2 = tmp;
end
factor = y;

end

function [VLong,VLat] = Rotate2(RotLong,RotLat,RotAng,VLong,VLat)
% '============================================================================
% '... Rotate subroutine      MTB  29 April 96
% ' This subroutine rotates a point through an angle around a pole, with
% ' the sense of rotation defined by the right-hand rule. Latitude is
% ' positive north (+Pi/2 --> -Pi/2) and longitude is positive east
% ' (+Pi --> -Pi or 0 --> 2Pi). Rotlat and rotlong represent latitude
% ' and longitude (or plunge and trend) of the rotation pole.
% ' Rotang is the rotation angle, vlat and vlong the latitude and
% ' longitude (or plunge and trend) of the vector to be rotated. The latitude
% ' and longitude (or plunge and trend) of the new vector are returned in
% ' radians as vlat and vlong. The program returns results consistent with the
% ' structural-geology convention that the sign of the plunge is positive when
% ' pointed down into the lower hemisphere (i.e. south hemisphere).
% '============================================================================
% DIM X1!(3), X!(3), A!(3, 3)
% '============================================================================
X(1) = cos(VLong) * cos(VLat);
X(2) = sin(VLong) * cos(VLat);
X(3) = sin(VLat);

C1 = cos(RotLong) * cos(RotLat);
C2 = sin(RotLong) * cos(RotLat);
C3 = sin(RotLat);

A1 = C1 * sin(RotAng / 2); 
U = C2 * sin(RotAng / 2);
V = C3 * sin(RotAng / 2); 
P = cos(RotAng / 2);

A(1, 1) = A1 * A1 - U * U - V * V + P * P;
A(2, 2) = U * U - V * V - A1 * A1 + P * P;
A(3, 3) = V * V - A1 * A1 - U * U + P * P;
A(2, 1) = 2 * (A1 * U + V * P); A(1, 2) = 2 * (A1 * U - V * P);
A(3, 1) = 2 * (V * A1 - U * P); A(1, 3) = 2 * (V * A1 + U * P);
A(3, 2) = 2 * (U * V + A1 * P); A(2, 3) = 2 * (U * V - A1 * P);
for I = 1:3
	X1(I) = A(I, 1) * X(1) + A(I, 2) * X(2) + A(I, 3) * X(3);
end

[VLong,VLat] = XYZ2TP(X1(1),X1(2),X1(3));

if abs(VLat) == pi/2
    VLong = 0;
end


end

function [betaI] = Betai(A,B,X)
% '=============================================================================
% '... Returns the incomplete beta function: Ix(a,b)
% '=============================================================================

if X < 0 || X > 1
    error('X not in range [0,1]')
end
if X == 0
    BT = 0;
else
    bt = exp(GammLn(A+B) - GammLn(A) - GammLn(B) + A * log(X) + B*log(1-X));
end
if X < (A+1)/(A+B+2)
    betaI = bt * Betacf(A,B,X)/A;
else
    betaI = 1 - bt*Betacf(B,A,1-X)/B;
end

if isnan(betaI)
    keyboard
end

end

function betaCF = Betacf(a,b,x)
% ... continued fraction fro Betai
itmax = 100;
eps = 0.0000003;
am = 1;
bm = 1;
az = 1;
qab = a+b;
qap = a+1;
qam = a-1;
bz = 1-qab*x/qap;
m = 0;
while m < itmax
    m = m+1;
    em = m;
    tem = 2*em;
    d = em*(b-m)*x/((qam+tem)*(a+tem));
    ap = az+d*am;
    bp = bz+d*bm;
    d = -(a+em)*(qab+em)*x/((a+tem)*(qap+tem));
    app = ap+d*az;
    bpp = bp+d*bz;
    aold = az;
    am = ap/bpp;
    bm = bp/bpp;
    az = app/bpp;
    bz = 1;
    if abs(az-aold) < eps*abs(az)
        m = 10*itmax;
%         return
    end
end
if abs(az-aold) >= eps*abs(az)
    error('a or b too bic or itmax too small');
else
    betaCF = az;
end

end

function y = GammLn(x)
% '=============================================================================
% '... Calculates the logrithm of the gamma function using the method of
% ' Lanczos. See Press and others (1986, p. 157) for details.
% '=============================================================================
y = gammaln(x);


end

function [trend,plunge] = xyz2tp(vec)
% function to return the trend and plunge of a 3D vector with components
% x,y,z. Vector need not be a unit vector, angles are returned in radians
% convention is x->y gives positive angle (e.g. NED reference frame)
% plunge is positive towards positive x (e.g. NED reference frame)

% trend [-pi,pi]
% plunge [-pi/2, pi/2]

x = vec(1);
y = vec(2);
z = vec(3);


trend = atan2(y,x);
mag = sqrt(x.^2+y.^2);
if mag == 0
    plunge = pi/2 * sign(z);
    trend = 0;
else
    plunge = atan(z./mag);
end

end

function [adjAngle] = adjang(angle,minAngle,range)
% '============================================================================
% '... The input Angle the variable to be converted. It can have any units.
% ' MinAng establishes the beginning of the desired angle range.
% ' For example, MinAng = -Pi/2 and Range = Pi would return a range:
% ' -Pi/2 to +Pi/2.  The newly adjusted angle is returned as adjAngle.
% Modified to MATLAB script by Chris Thissen, 2014
% '============================================================================
adjAngle = angle - range*floor((angle-minAngle)/range);

end

function [VLong,VLat] = rotate_yale(RotLong,RotLat,RotAng,VLong,VLat)
% rotate a vector given by [VLong,VLat] about an axis defined
% by [RotLong,RotLat] by and angle RotAng. Makes use of the rodriguez
% rotation formula

%... construct unit vectors
[rotX,rotY,rotZ] = sph2cart(RotLong,RotLat,1);
[vX,vY,vZ] = sph2cart(VLong,VLat,1);
v = [vX,vY,vZ];
k = [rotX,rotY,rotZ];

vRot = v*cos(RotAng) + cross(k,v)*sin(RotAng) + k*dot(k,v)*(1-cos(RotAng));

[VLong,VLat,~] = cart2sph(vRot(1),vRot(2),vRot(3));
end




