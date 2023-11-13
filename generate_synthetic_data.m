function [x,ttX,ttY,LagVector,Coeff,Delta]=generate_synthetic_data(n,t,varargin)
%----------------------------------------------------------------------maca
%------------------------------------------------------------------<header>
% This function generates synthetic data for analysis.
% Inputs:
%   n: Number of data points to generate
%   t: Time parameter affecting the data generation
%   varargin: Additional optional parameters
%-----------------------------------------------------------------<\header>

% Initialize default values for various parameters and set up an input parser
% These defaults control various aspects of the synthetic data generation process

defaultSet = [0;0;0];
defaultFull = repmat([0;0;0],1,t);
defaultDeltaLower=zeros(1,9);
defaultDeltaUpper=50*ones(1,9);
defaultCoeffSigma = 0;
defaultXnoise=[1,1,1];
defaultYnoise=0.05;
defaultLag=5;
% Create an input parser object to handle and validate input arguments
P=inputParser;
% Define default values and/or validators for various parameters
validScalarNonnegNum = @(x) isnumeric(x) && all(all(x >= 0));
validScalarPosNum = @(x) isnumeric(x) && all(all(x > 0));
validFullSet= @(x) validateCoeffcients(x) && all(all(size(x)==[3,t]));
addParameter(P,'fullSet',   defaultFull,validFullSet)
addParameter(P,'coeffSigma',defaultCoeffSigma,validScalarNonnegNum)
addParameter(P,'deltaUpper',defaultDeltaUpper,@(x) isnumeric(x))
addParameter(P,'deltaLower',defaultDeltaLower,@(x) isnumeric(x))
addParameter(P,'Xnoise',    defaultXnoise,validScalarNonnegNum)
addParameter(P,'Ynoise',    defaultYnoise,validScalarNonnegNum)
addParameter(P,'Outliers',  false)
addParameter(P,'startSet',  defaultSet,validScalarNonnegNum)
addParameter(P,'Verbose',   false)
addParameter(P,'lagNum',    defaultLag, validScalarPosNum)
% Parse the input arguments
parse(P,varargin{:});
% Extract the values from parsed input
c0=P.Results.startSet;
if length(c0)==1
    c0=[c0;c0;c0];
end
fullSet=P.Results.fullSet;
Verbose=P.Results.Verbose;
Outliers=P.Results.Outliers;
NoiseY=P.Results.Ynoise;
NoiseX=P.Results.Xnoise;
coeffSigma=P.Results.coeffSigma;
Delta.Lower=P.Results.deltaLower;
Delta.Upper=P.Results.deltaUpper;
lag=P.Results.lagNum;
% Initialize parameters for the data generation process
expressedVariables=9;
LagVector=fliplr(lag*(1:expressedVariables));
lag_max_lead=lag+lag*expressedVariables;
Delta.max_lag=lag_max_lead;
Delta.Target=LagVector;
t_safety=101*t;
safety_tail=max(101*t,lag_max_lead+1);
k=n+t_safety+safety_tail;
coeff_orig=[0.2028   -0.8473    0.4909];
coeff_orig=coeff_orig./norm(coeff_orig);
x1_coeff=circshift(coeff_orig,0); %[ 0.2028  -0.8473   0.4909];
x2_coeff=circshift(coeff_orig,1); %[ 0.4909   0.2028  -0.8473];
x3_coeff=circshift(coeff_orig,2); %[-0.8473   0.4909   0.2028];
lin_coeff=[-0.7071,0.7071] + mvnrnd(zeros(1,2),coeffSigma*eye(2));
lin_coeff=lin_coeff./norm(lin_coeff);%Norm==1
int_coeff=[-0.2300    0.2300   -0.7146   -0.4732    0.3995]+ mvnrnd(zeros(1,5),coeffSigma*diag([0.1 0.1 1 1 1]));
int_coeff=int_coeff./norm(int_coeff);
% safe_n=n+safety_tail;
% kernel_coeff=((1:safe_n).^2.*sin((1:safe_n).^2))./(1+(1:safe_n).*cos((1:safe_n)));
% kernel_coeff=kernel_coeff+ mvnrnd(zeros(1,size(kernel_coeff,2)),std(kernel_coeff)*coeffSigma*eye(size(kernel_coeff,2)))/norm(kernel_coeff);
% sigma_kernel=37/23;
% gamma_kernel=1/(2*(sigma_kernel)^2);
% Define different models for data generation
linear=@(X)normalize(X(:,1)*lin_coeff(1)...
    +X(:,2)*lin_coeff(2),'medianiqr');
interactions=@(X)normalize(X(:,1)*int_coeff(1)...
    +X(:,2)*int_coeff(2)...
    +X(:,1).^2*int_coeff(3)...
    +X(:,2).^2*int_coeff(4)...
    +X(:,1).*X(:,2)*int_coeff(5),'medianiqr');
% nonLinear=@(X)normalize(sign((X(:,1))*lin_coeff(1)).*log2(abs(X(:,1)*lin_coeff(1)))...
%     +sign((X(:,2))*lin_coeff(2)).*nthroot(abs(X(:,2)*lin_coeff(2)),3),'medianiqr');
% kernel=@(X)normalize(X*kernel_coeff','medianiqr');
% exponential=@(X)normalize(exp(nthroot(X(:,1)*lin_coeff(1),5).^2 ...
%     +nthroot(X(:,2)*lin_coeff(2),5).^2));
powerRatio=@(X)normalize(power(((X(:,2))*lin_coeff(2)),2)...
    ./(1+power(((X(:,1)).*lin_coeff(1)),3)),'medianiqr');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filter_stable = false;
c1=zeros(t,1);
c2=zeros(t,1);
c3=zeros(t,1);
% Check if fullSet is the default set, otherwhise generate stable AR
% coefficients or check stablity of the provided coeff
if all(fullSet==defaultSet)
    count=0;
    while ~filter_stable
        count=count+1;
        if all(c0==defaultSet)
            c0=rand(3,1);
        end
        c1(1)=c0(1);
        c2(1)=c0(2);
        c3(1)=c0(3);
        c1(2:t)=randperm(100,t-1)*0.001*(1-c1(1));
        c2(2:t)=randperm(100,t-1)*0.001*(1-c1(1));
        c3(2:t)=randperm(100,t-1)*0.001*(1-c1(1));
        coefficients1 = [1; -c1];
        coefficients2 = [1; -c2];
        coefficients3 = [1; -c3];
        if max(abs(roots(coefficients1))) < 1 && max(abs(roots(coefficients2))) < 1 && max(abs(roots(coefficients3))) < 1
            filter_stable = true;
            a1=c1;
            a2=c2;
            a3=c3;
        end
        if count >n*t*1000
            msg="No stable filter configurations."...
                +newline+"Try:"...
                +newline+" - Different starting points"...
                +newline+" - Different seed"...
                +newline+" - Lower autoregression order"...
                ;
            error(msg)
        end
    end
else
    a1=fullSet(1,:);
    a2=fullSet(2,:);
    a3=fullSet(3,:);
end
Coeff=[a1';a2';a3'];
if Verbose
    CoeffStr="";
    for i=1:3
        CoeffStr=CoeffStr+string(i)+": ";
        for j=1:t
            CoeffStr=CoeffStr+string(Coeff(i,j))+", ";
        end
        CoeffStr=CoeffStr+newline;
    end
    msg="Generating data with AR order "+string(t)...
        + newline +CoeffStr;
    disp(msg);
end
% Generate base data using multivariate normal distribution; force zero
% correlation between the three base variables
base=mvnrnd(zeros(3,1),eye(3),k);
% Apply filters to base data to create X1, X2, X3
X1=filter(1,[1,-a1],base(:,1));
X2=filter(1,[1,-a2],base(:,2));
X3=filter(1,[1,-a3],base(:,3));
X1=X1(t_safety+1:end,:);
X2=X2(t_safety+1:end,:);
X3=X3(t_safety+1:end,:);
X1=normalize(X1,'medianiqr');
X2=normalize(X2,'medianiqr');
X3=normalize(X3,'medianiqr');
x=zeros(size(X2,1),9);
X1_x=X1*x1_coeff;
X2_x=X2*x2_coeff;
X3_x=X3*x3_coeff;
if Verbose
    disp(char(hex2dec('2713')))
end
% Create matrix X
for i=1:3
    XNoise=mvnrnd(zeros(3,1),diag(NoiseX),size(x,1));
    x1=(i*3-2);
    x2=(i*3-1);
    x3=(i*3);
    %"Measurement" error on the measured x
    x(:,x1)=X1_x(:,i)+XNoise(:,1);
    x(:,x2)=X2_x(:,i)+XNoise(:,2);
    x(:,x3)=X3_x(:,i)+XNoise(:,3);
end
% Insert outliers into the data if required
if Outliers
    if Verbose
        msg="Adding Outliers...";
        disp(msg)
    end
    OUTProb1=0.05;
    OUTProb2=[3/6,2/6,1/6];
    OuliersCount=0;
    span_outliers=lag_max_lead+1:size(x,1)-(safety_tail-lag_max_lead);
    for i=span_outliers
        r=binornd(1,OUTProb1);
        if r
            OuliersCount=OuliersCount+1;
            r2=(1:3).*mnrnd(1,OUTProb2,1);
            r2(r2==0)=[];
            if  r2==3
                P=randperm(9,3);
                x(i,P)=x(i,P)+sign(x(i,P)-median(x(:,P)))*1.*iqr(x(:,P));
                
            elseif r2==2
                P=randperm(9,2);
                x(i,P)=x(i,P)+sign(x(i,P)-median(x(:,P)))*2.*iqr(x(:,P));
                
            elseif r2==1
                P=randperm(9,1);
                x(i,P)=x(i,P)+sign(x(i,P)-median(x(:,P)))*3*iqr(x(:,P));
                
            end
        end
    end
    totalOutliers=(OuliersCount/n)*100;
    if Verbose
        
        msg=string(sprintf('%2.1f',round(totalOutliers,1)))+"% of total observations";
        disp(msg)
        disp(char(hex2dec('2713')))
    end
end
% Insert lag into the data
j=1;
for i=fliplr(1:9)
    VariableLag=-lag*j;
    x(:,i)=circshift(x(:,i),VariableLag);
    j=j+1;
end

% Create Y values using different models
D=([X1,X2]);
yNoise=NoiseY*mvnrnd(zeros(6,1),eye(6),size(D,1));
% DD=squareform(pdist(D,'squaredeuclidean'));
% gauss_kernel=exp(-gamma_kernel*sqrt(DD));
y_lin=linear(D)+yNoise(:,1);
y_int=interactions(D)+yNoise(:,2);
% y_nlin=nonLinear(D)+yNoise(:,3);
% y_kernel=kernel(gauss_kernel)+yNoise(:,4);
% y_exp=exponential(D)+yNoise(:,5);
y_enzyme= powerRatio(D)+yNoise(:,6);
% Create timetables for X and Y data
time_idx=datetime('yesterday')+minutes(0:lag:lag*(size(D,1)-1));
ttX = timetable('Size',[size(D,1),9],'VariableTypes',repmat("double",[1,9]),'RowTimes',time_idx);
ttY = timetable('Size',[size(D,1),3],'VariableTypes',repmat("double",[1,3]),'RowTimes',time_idx);
ttX.Properties.VariableNames={'X1_1','X2_1','X3_1','X1_2','X2_2','X3_2','X1_3','X2_3','X3_3'};
% ttY.Properties.VariableNames={'y_lin','y_interactions','y_kernel','y_nlin','y_exp','y_enzyme'};
ttY.Properties.VariableNames={'y_lin','y_interactions','y_enzyme'};
% Adjust the lag and augument the matrix
ttX{:,:}=x;
ttY{:,:}=[y_lin,y_int,y_enzyme];
[ttX,ttY]=timetable_augument(ttX,ttY,Delta);
if  max(Delta.Upper)>lag_max_lead
    lag_max_lead=max(Delta.Upper);
end
if  min(Delta.Lower)<0
    lag_max_lead=lag_max_lead-min(Delta.Lower);
end
Delta.max_lag=lag_max_lead;
% Trim and finalize the data
ttY=ttY(1:end-(safety_tail-lag_max_lead),:);
ttX=ttX(1:end-(safety_tail-lag_max_lead),:);
if Verbose
    msg="done";
    disp(msg)
    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%HELPER FUNCTIONS%%
% Additional function to validate coefficients
function flag=validateCoeffcients(x)
c1=x(1,:)';
c2=x(2,:)';
c3=x(3,:)';
coefficients1 = [1; -c1];
coefficients2 = [1; -c2];
coefficients3 = [1; -c3];
if max(abs(roots(coefficients1))) < 1 && max(abs(roots(coefficients2))) < 1 && max(abs(roots(coefficients3))) < 1
    flag = true;
else
    flag= false;
    warning('coefficients not stable')
end

end
