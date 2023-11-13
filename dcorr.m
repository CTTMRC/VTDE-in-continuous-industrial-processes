function dCorr=dcorr(x,y)
%----------------------------------------------------------------------maca
%------------------------------------------------------------------<header>
% Function to calculate the distance correlation between two data sets x and y.
    % Inputs:
    %   x - First dataset (matrix of observations)
    %   y - Second dataset (matrix of observations)
    % Output:
    %   dCorr - distance correlation
% Distance correlation algorithm published in: 
% https://doi.org/10.1214/009053607000000505
%-----------------------------------------------------------------<\header>

% Determine the number of observations
n=size(x,1);
linIdx=sub2ind([n,n],1:n,1:n);
% Calculate pairwise distances for dataset y and mean-adjust them
SF_Y=squareform(pdist(y));
Bbar=sum(sum(SF_Y))/(n-1)/(n-2);
Bk=repmat(sum(SF_Y,2)/(n-2),1,size(SF_Y,2));
Bl=repmat(sum(SF_Y,1)/(n-2),size(SF_Y,1),1);
EY=Bk+Bl-Bbar;
EY(linIdx)=0;
Bkl=SF_Y-EY;
% Calculate pairwise distances for dataset x and mean-adjust them
SF_X=squareform(pdist(x));
Abar=sum(sum(SF_X))/(n-1)/(n-2);
Ak=repmat(sum(SF_X,2)/(n-2),1,size(SF_X,2));
Al=repmat(sum(SF_X,1)/(n-2),size(SF_X,1),1);
EX=Ak+Al-Abar;
EX(linIdx)=0;
Akl=SF_X-EX;
% Calculate the components of distance covariance and variance
V2XY=sum(sum(Akl.*Bkl'));
V2X=sum(sum(Akl.*Akl'));
V2Y=sum(sum(Bkl.*Bkl'));
% Calculate the distance correlation
if V2X>0&&V2Y>0
    sqDCorr=abs(V2XY)./real(sqrt(V2X*V2Y'));
else
    sqDCorr=0;
end
% Final distance correlation value
dCorr=real(sqrt(sqDCorr));
end
