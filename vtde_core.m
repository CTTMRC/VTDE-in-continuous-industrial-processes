reps=3;
ffFactors=fullfact([3,3,2,2]);
numFactors=length(ffFactors);
sizeFactors=cellstr(["100___","1000__","10000_"])';
complexityFactors=cellstr(["linear______","interactions","power_ratio_"])';
arOrder=cellstr(["_1","_2"])';
arStrenght= cellstr(["_high","_low"])';
rowNames=strcat(sizeFactors(ffFactors(:,1)),complexityFactors(ffFactors(:,2)),arOrder(ffFactors(:,4)),arStrenght(ffFactors(:,3)));
methods=cellstr("DistanceCorrelation");
repname="";
for z=1:reps
    repname=([repname,"_rep"+num2str(z)]);
end
repname=repname(2:end);
ResultsTruthfulness=table('Size',[numFactors*reps,numel(methods)],'VariableTypes',repmat("doublenan",[1,numel(methods)]),'VariableNames',methods,'RowNames',repmat(rowNames,[reps,1])+repelem(repname',numFactors));
sizeField=cellstr(["hundred","thousand","tenThousand"]);
complexityField=cellstr(["C1","C2","C3"]);
arField1=(["H","L"]);
arField2=(["1","2"]);
tic
vtdeDCorrelation=   zeros(numFactors*reps,9);
timeDCorrelation=   zeros(numFactors*reps,1);
%%
t1=tic;
progress=0;
wait_bar=waitbar(0,[num2str(progress) , '%']);
for repetition=1:reps
    [ExperimentData,LagVect,Delta]=generate_experiment_set('coeffSigma',0.01,'rngSeed',repetition);
    buffer=numFactors*(repetition-1);
    for experiment=12%:numFactors
        ttXCalib=   ExperimentData.(arField1(ffFactors(experiment,3))+arField2(ffFactors(experiment,4)))...
            .(sizeField{ffFactors(experiment,1)}).calib.TT;
        yCalib=     ExperimentData.(arField1(ffFactors(experiment,3))+arField2(ffFactors(experiment,4)))...
            .(sizeField{ffFactors(experiment,1)}).calib...
            .(complexityField{ffFactors(experiment,2)});
        ttXCalib=ttXCalib(:,(std(ttXCalib{:,:},1))>0);
        tic
        [vtdeDCorrelation(experiment+buffer,:),distanceDCORR]=vtde_compute_dcorr(ttXCalib,yCalib,Delta);
        timeDCorrelation(experiment+buffer)=toc;
        ResultsTruthfulness{experiment+buffer,"DistanceCorrelation"}=distanceDCORR;
        progress=progress+1;
        waitbar(progress/(numFactors*reps),wait_bar,[num2str(progress) , '%']);
    end
end
totalTime=toc(t1);
%% Analysis
SS = categorical(repmat(ffFactors(:,1),numel(methods)*reps,1),[1 2 3],{'100' '1000' '10000'});
CPLX = categorical(repmat(ffFactors(:,2),numel(methods)*reps,1),[1 2 3],{'linear' 'interaction' 'power'});
AROrder = categorical(repmat(ffFactors(:,4),numel(methods)*reps,1),[1 2 ],cellstr(["1","2"]));
ARStrength = categorical(repmat(ffFactors(:,3),numel(methods)*reps,1),[1 2 ],cellstr(["H","L"]));
groups= {(repelem(ResultsTruthfulness.Properties.VariableNames,numFactors*reps)');...
    SS;...
    CPLX;...
    AROrder;...
    ARStrength}...
    ;
[truth,~,statsTruth] = anovan(reshape(ResultsTruthfulness{:,:},[],1),...
    groups,'model','interaction','varnames',{'method','size','complexity','ARStrength','AROrder'});




