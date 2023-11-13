function [vtde,euclideanDistance]=ffAnalysisMRmr(X,Y,Delta,varargin)
defaultOptimOptions=struct(...
    'nPop',500,...
    'maxGeneration',50,...
    'CR',0.01,...
    'MR',0.02,...
    'PRESSURE',4,...
    'Elitism',0.01...
    );
P=inputParser;
addRequired(P,'Delta')
addOptional(P,'Options',defaultOptimOptions,@(x)validOptionStructure(x,defaultOptimOptions));
parse(P,Delta);
upper=Delta.Upper;
lower=Delta.Lower;
target=Delta.Target;
maxLag=upper-lower;
span=[0,cumsum(maxLag)];
spanReal=span(1:length(upper));
optimOptions=P.Results.Options;
nPop=optimOptions.nPop;
maxGeneration=optimOptions.maxGeneration;
freshStartRate=optimOptions.CR;
mutationRate=optimOptions.MR;
mutationPressure=optimOptions.PRESSURE;
elitism=optimOptions.Elitism;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask=~isoutlier(Y,'quartiles');
X=X{mask,:};
Y=Y(mask,:);
S_max=-2e20;
S_min=2e20;
[~,~,MIxy,MIxx]=statistics4(X',Y');
% [L,~]=size(MIxx);
C=length(maxLag);
redundancyWeight=1;
%starting genes:
population= floor(maxLag'.*rand(C,nPop))+1;
Sh=nan(1,nPop);
Sh2=nan(1,nPop);
Sh(1)=0;
Sh2(1)=0.05;
%**************************************************************************
%Loop over the generations:
%**************************************************************************
generation=0;
while generation<=maxGeneration&&abs(Sh(generation+1)-Sh2(generation+1))>0.0002
    generation=generation+1;
    %ranking:
    S=zeros(1,nPop);
    for k=1:nPop
        Xind=population(:,k)+spanReal';
        %relevance:
        mean_MIxy=mean(sum(MIxy(Xind)));
        %redundancy:
        mean_MIxx=0;
        for n=1:C
            for m=1:C
                mean_MIxx=mean_MIxx+MIxx(Xind(n),Xind(m));
            end
        end
        mean_MIxx=(mean_MIxx-C)/(C^2-C);
        %fitness:
        S(k)=mean_MIxy-redundancyWeight*mean_MIxx;
        if S(k)>S_max
            S_max=S(k);
        end
        if S(k)<S_min
            S_min=S(k);
        end
    end
    Sh(generation+1)=mean(S);
    Sh2(generation+1)=max(S);
    [~,sortIndex]=sort(S,'descend');
    rearrangedPopulation=population(:,sortIndex);
  
    if generation>1
        eliteIndividuals=rearrangedPopulation(:,1:ceil(elitism*nPop));
    end
    %crossover:
    for k=1:nPop
        r=rand(2,1);
        father=round((nPop-1)*(exp(mutationPressure*r(1))-1)/(exp(mutationPressure)-1))+1;
        mother=round((nPop-1)*(exp(mutationPressure*r(2))-1)/(exp(mutationPressure)-1))+1;
        r=normrnd(0,1,1,9);
        population(r<=0,k)=rearrangedPopulation(r<=0,mother);
        population(r>0,k)=rearrangedPopulation(r>0,father);
        %mutation
        mutationChance=rand;
        if mutationChance<mutationRate
            mutationVector=rand(C,1) ;
            [~,mutationLow]=min(mutationVector);
            [~,mutationHigh]=max(mutationVector);
            mutatedPop=rearrangedPopulation(:,k);
            mutatedPop(mutationLow,:)=rearrangedPopulation(mutationHigh,k);
            mutatedPop(mutationHigh,:)=rearrangedPopulation(mutationLow,k);
            population(:,k)=mutatedPop;
        end
        %Verifying not used features:
        unusedFeatures=false(C,max(maxLag));
        for n=1:max(maxLag)
            for i=1:C
                usedFeatureCheck=find(population(i,:)==n);
                features=[usedFeatureCheck,0];
                
                if sum(features)==0
                    unusedFeatures(i,n)=true;
                end
            end
        end
        nUnused=zeros(1,C);
        if sum(sum(unusedFeatures))>0
            for i=1:C
                unused=find(unusedFeatures(i,:));
                if sum(unused)>0
                    nUnused(i)=sum(unusedFeatures(i,:));
                    population(i,end+1-nUnused(i):end)=unused;
                end
            end
        end
    end
    if generation>1
        eliteIndex=sortIndex(end+1-max(max(nUnused),0)-ceil(elitism*nPop):end-max(max(nUnused),0));
        population(:,eliteIndex)=eliteIndividuals;
    end
    freshStartChance=rand;
    if freshStartChance<freshStartRate
        freshStartIndex=sortIndex(end-(max(max(nUnused),0)+1+ceil(elitism*nPop)));
        population(:,freshStartIndex)=floor(maxLag'.*rand(C,1))+1;
    end
end
vtde=rearrangedPopulation(:,1)-1;
distance=(vtde-(target-lower))./maxLag;
relevantDistance=(distance([1:2,4:5,7:8]));
euclideanDistance=norm([relevantDistance,zeros(size(relevantDistance))]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

function OK=validOptionStructure(in,valid)
valid=fields(valid);
if isstruct(in)
    idx=zeros(1,length(in));
    check=fields(in);
    for i=1:length(in)
        idx(i)=contains(check{i},valid);
    end
    OK=all(idx);
else
    OK=false;
end
end