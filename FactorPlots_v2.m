close all
load C:\Users\marco.cattaldo\MATLAB\Projects\ffAnalysisVTDE\Data\RESULTS20230120.mat
ffFactors=fullfact([3,3,4]);
numFactors=size(ffFactors,1);
methodNames=[{'r'};
    {'\rho'};
    {'\tau'};
    {'PLS-SR'};
    {'PLS-Beta'};
    {'VTR-LS'};
    {'MIC'};
    {'N-Norm'};
    {'MI'};
    {'dCorr'};
    {'GA-MI'};
    {'KPLS'};
    {'RF'};
    {'VTR-GMM'}];

methodCategories=cellstr(["MIC"
    "dCorr"
    "GA-MI"
    "N-Norm"
    "MI"
    "RF"
    "KPLS"
    "\rho"
    "\tau"
    "PLS-SR"
    "PLS-Beta"
    "r"
    "VTR-LS"
    "VTR-GMM"]);
%ResultsTruthfulness.Properties.VariableNames
groups= {(repelem(methodNames,numFactors*3));...
    repmat(repmat(ffFactors(:,1),14,1),3,1);...
    repmat(repmat(ffFactors(:,2),14,1),3,1);...
    repmat(repmat(ffFactors(:,3),14,1),3,1)}...
    ;
%%
SS = categorical(repmat(ffFactors(:,1),3,1),[1 2 3],{'100' '1000' '10000'});
CPLX = categorical(repmat(ffFactors(:,2),3,1),[1 2 3],{'linear' 'interaction' 'power'});
AR = categorical(repmat(ffFactors(:,3),3,1),[1 2 3 4],cellstr(["H1","L1","H2","L2"]));
ARNum = categorical(repmat(ffFactors(:,3),3,1),[1 2 3 4],cellstr(["1","1","2","2"]));
ARLvl = categorical(repmat(ffFactors(:,3),3,1),[1 2 3 4],cellstr(["H","L","H","L"]));
MethodGroup = nan(length(groups{1}),1);
MethodGroup(contains(groups{1},'DistanceCorrelation') |...
    contains(groups{1},'MI') |...
    contains(groups{1},'MIC') |...
    contains(groups{1},'GA-MI') |...
    contains(groups{1},'N-Norm'))=1;
MethodGroup(contains(groups{1},'PLS-SR')|...
    contains(groups{1},'PLS-Beta') |...
    contains(groups{1},'r') |...
    contains(groups{1},'\tau') |...
    contains(groups{1},'\rho') |...
    contains(groups{1},'VTR-LS')) =2;
MethodGroup(contains(groups{1},'RF') |...
    contains(groups{1},'KPLS'))=3;
MethodGroup(contains(groups{1},'VTR-GMM'))=4;
%%
XX = ResultsTruthfulness{:,:};
TT = table(reshape(XX,[],1));
TT.Properties.VariableNames(1)="XX";
TT.Method = groups{1};
TT.Method=replace(TT.Method,"_","-");
TT.Method=replace(TT.Method,"MI-Index","MI");
TT.Method=replace(TT.Method,"GA-MRmr","GA-MI");
TT.Size = repmat(SS,14,1);
TT.Complexity   = repmat(CPLX,14,1);
TT.ARnum =  repmat(ARNum,14,1);
TT.ARlvl = repmat(ARLvl,14,1);
TT.MethodGroup = MethodGroup;
TT.Hits=TT.XX==0;
TT.Time=allTIME;
[P,T,STATS] = anovan(TT.XX,{TT.Size TT.Complexity TT.ARnum TT.ARlvl TT.Method},"model",3,'varnames',{'Size','Complexity','DynOrder', 'DynStr' ,'Method'});
Tssq=cellfun(@(x)x/T{end,2},T(2:end,2));
T(1,8)={'Explained  Variance [%]'};
T(2:end,8)=num2cell(Tssq*100);
emptyMask=cellfun(@(x)isempty(x),T(2:end,2:end));
T(28,5:7)=[{0},{0},{0}];
T(27,6:7)=[{0},{0}];
T=cell2table(T(2:end,:),'VariableNames',T(1,:));
T{:,2:end}=round(T{:,2:end},2);
T.("Singular?")=[];
%%
figure();
[~,~,f]=multcompare(STATS,'Dimension',5,'CType','lsd');
plotData=readMultcomparePlot(f);
[~,I]=sortrows(plotData,"m",'ascend');
performanceOrder=double(categorical(plotData.Method(I)));
uniqueMethods=unique(plotData.Method);
catMethods=categorical(plotData.Method,methodCategories);
f3=figure();
T1=tiledlayout(1,1);
tmp = (plotData(:,:));%,:);%
t5=nexttile;
h=scatter(double(catMethods),tmp.ul); legend off
h.XData=h.XData(:,I);
h.YData=h.YData(:,I);
h.Marker = '_';
h.SizeData= 72;
h.MarkerFaceColor = "r";
h.MarkerEdgeColor = "r";
UL=h.YData;
hold on
h=scatter(double(catMethods),tmp.ll); legend off
h.XData=h.XData(:,I);
h.YData=h.YData(:,I);
h.Marker = '_';
h.SizeData= 72;
h.MarkerFaceColor = "r";
h.MarkerEdgeColor = "r";
for j=1:size(h(1).XData,2)
    line([h.XData(j),h.XData(j)],[UL(j),h.YData(j)],'Color','k')
end
h=scatter(double(catMethods),tmp.m); legend off
h.XData=h.XData(:,I);
h.YData=h.YData(:,I);
h.Marker = "o";
h.SizeData= 54;
h.MarkerFaceColor = "r";
h.MarkerEdgeColor = 'k';
t5.Title.String='General method perfomance';
t5.Title.FontSize=26;
t5.YLim=[0 1];
t5.XLim=[0.5 14.5];
t5.XTick=1:14;
t5.XTickLabel=categories(catMethods);
t5.XTickLabelRotation=30;
t5.PlotBoxAspectRatio=[3,1,1];
t5.FontSize=26;
box on;
T1.XLabel.String="Method";T1.XLabel.FontWeight="bold";T1.XLabel.FontSize=26;
T1.YLabel.String="D_{method}"+newline;T1.YLabel.FontWeight="bold";T1.YLabel.FontSize=26;

%%
figure();
[~,~,f]=multcompare(STATS,'Dimension',[1,2,5],'CType','lsd');
plotData=readMultcomparePlot(f);
S=cellstr(["v","s","o","x","^","d"]);
C=cellstr(["r","b","#017d03","k"]);
MSz=8;
f1=figure();
T1=tiledlayout(3,1);
T1.TileSpacing = 'tight';
T1.Padding = 'compact';
t1=nexttile(1);
tmp = (plotData(plotData.Size=="100",:));%,:);%
catMethods=categorical(tmp.Method,methodCategories);
offset=[-0.2,0,0.2];
UL=nan(3,14);
h=gscatter(double(catMethods),tmp.ll,tmp.Complexity); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    UL(i,:)=h(i).YData;
end
hold on
h=gscatter(double(catMethods),tmp.ul,tmp.Complexity); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    for j=1:size(h(1).XData,2)
        line([h(i).XData(j),h(i).XData(j)],[UL(i,j),h(i).YData(j)],'Color','k')
    end
end
h=gscatter(double(catMethods),tmp.m,tmp.Complexity); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = S(i);
    h(i).MarkerSize= MSz;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = 'k';%C(i) ;
    h(i).Color = 'k';%C(i);
end
t1.Title.String='Size = 100';
t1.Title.FontSize=26;
t1.YLim=[0 1.1];
t1.XLim=[0.5 14.5];
t1.XTick=1:14;
t1.XTickLabel=categories(catMethods);
t1.XTickLabelRotation=30;
t1.FontSize=18;
box on;
t3=nexttile(2);
tmp = (plotData(plotData.Size=="1000",:));%,:);%
catMethods=categorical(tmp.Method,methodCategories);
offset=[-0.2,0,0.2];
UL=nan(3,14);
h=gscatter(double(catMethods),tmp.ll,tmp.Complexity); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    UL(i,:)=h(i).YData;
end
hold on
h=gscatter(double(catMethods),tmp.ul,tmp.Complexity); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    for j=1:size(h(1).XData,2)
        line([h(i).XData(j),h(i).XData(j)],[UL(i,j),h(i).YData(j)],'Color','k')
    end
end
h=gscatter(double(catMethods),tmp.m,tmp.Complexity); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = S(i);
    h(i).MarkerSize= MSz;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = 'k';%C(i) ;
    h(i).Color = 'k';%C(i);
end
t3.Title.String='Size = 1000';
t3.Title.FontSize=26;
t3.YLim=[0 1.1];
t3.XLim=[0.5 14.5];
t3.XTick=1:14;
t3.XTickLabel=categories(catMethods);
t3.XTickLabelRotation=30;
t3.FontSize=18;
box on;
t5=nexttile(3);
tmp = (plotData(plotData.Size=="10000",:));%,:);%
catMethods=categorical(tmp.Method,methodCategories);
offset=[-0.2,0,0.2];
UL=nan(3,14);
h=gscatter(double(catMethods),tmp.ll,tmp.Complexity); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    UL(i,:)=h(i).YData;
end
hold on
h=gscatter(double(catMethods),tmp.ul,tmp.Complexity); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    for j=1:size(h(1).XData,2)
        line([h(i).XData(j),h(i).XData(j)],[UL(i,j),h(i).YData(j)],'Color','k')
    end
end
h=gscatter(double(catMethods),tmp.m,tmp.Complexity); legend off
for i=1:length(h)
    h(i).XData=h(i).XData(:,performanceOrder);
    h(i).YData=h(i).YData(:,performanceOrder);
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = S(i);
    h(i).MarkerSize= MSz;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = 'k';%C(i) ;
    h(i).Color = 'k';%C(i);
end
t5.Title.String='Size = 10000';
t5.Title.FontSize=26;
t5.YLim=[0 1.1];
t5.XLim=[0.5 14.5];
t5.XTickLabelRotation=30;
t5.XTickLabel=categories(catMethods);
t5.FontSize=26;
box on;
leg=legend(h,'Orientation', 'vertical');
leg.Layout.Tile = 'east';
leg.Interpreter = 'tex';
leg.Title.String="Complexity";
box on;
T1.XLabel.String="Method";T1.XLabel.FontWeight="bold";T1.XLabel.FontSize=26;
T1.YLabel.String="D_{method}"+newline;T1.YLabel.FontWeight="bold";T1.YLabel.FontSize=26;

%%
figure();
[~,~,f]=multcompare(STATS,'Dimension',[1,2,5],'CType','lsd');
plotData=readMultcomparePlot(f);
S=cellstr(["v","s","o","x","^","d"]);
C=cellstr(["r","b","#017d03","k"]);
MSz=8;
f2=figure();
T1=tiledlayout(3,1);
T1.TileSpacing = 'tight';
T1.Padding = 'compact';
t1=nexttile(1);
tmp = (plotData(plotData.Complexity=="linear",:));%,:);%
catMethods=categorical(tmp.Method,methodCategories);
offset=[-0.2,0,0.2];
UL=nan(3,14);
h=gscatter(double(catMethods),tmp.ll,tmp.Size); legend off
for i=1:length(h)
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    UL(i,:)=h(i).YData;
end
hold on
h=gscatter(double(catMethods),tmp.ul,tmp.Size); legend off
for i=1:length(h)
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    for j=1:size(h(1).XData,2)
        line([h(i).XData(j),h(i).XData(j)],[UL(i,j),h(i).YData(j)],'Color','k')
    end
end
h=gscatter(double(catMethods),tmp.m,tmp.Size); legend off
for i=1:length(h)
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = S(i);
    h(i).MarkerSize= MSz;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = 'k';%C(i) ;
    h(i).Color = 'k';%C(i);
end
t1.Title.String='Complexity = linear';
t1.Title.FontSize=26;
t1.YLim=[0 1.1];
t1.XLim=[0.5 14.5];
t1.XTick=1:14;
t1.XTickLabel=categories(catMethods);
t1.XTickLabelRotation=30;
t1.FontSize=26;
box on;
t3=nexttile(2);
tmp = (plotData(plotData.Complexity=="interaction",:));%,:);%
catMethods=categorical(tmp.Method,methodCategories);
offset=[-0.2,0,0.2];
UL=nan(3,14);
h=gscatter(double(catMethods),tmp.ll,tmp.Size); legend off
for i=1:length(h)
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    UL(i,:)=h(i).YData;
end
hold on
h=gscatter(double(catMethods),tmp.ul,tmp.Size); legend off
for i=1:length(h)
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    for j=1:size(h(1).XData,2)
        line([h(i).XData(j),h(i).XData(j)],[UL(i,j),h(i).YData(j)],'Color','k')
    end
end
h=gscatter(double(catMethods),tmp.m,tmp.Size); legend off
for i=1:length(h)
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = S(i);
    h(i).MarkerSize= MSz;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = 'k';%C(i) ;
    h(i).Color = 'k';%C(i);
end
t3.Title.String='Complexity = interactions';
t3.Title.FontSize=26;
t3.YLim=[0 1.1];
t3.XLim=[0.5 14.5];
t3.XTick=1:14;
t3.XTickLabel=categories(catMethods);
t3.XTickLabelRotation=30;
t3.FontSize=26;
box on;
t5=nexttile(3);
tmp = (plotData(plotData.Complexity=="power",:));%,:);%
catMethods=categorical(tmp.Method,methodCategories);
offset=[-0.2,0,0.2];
UL=nan(3,14);
h=gscatter(double(catMethods),tmp.ll,tmp.Size); legend off
for i=1:length(h)
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    UL(i,:)=h(i).YData;
end
hold on
h=gscatter(double(catMethods),tmp.ul,tmp.Size); legend off
for i=1:length(h)
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = '_';
    h(i).MarkerSize= MSz*2;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = C{i} ;
    h(i).Color = 'k';%C(i);
    for j=1:size(h(1).XData,2)
        line([h(i).XData(j),h(i).XData(j)],[UL(i,j),h(i).YData(j)],'Color','k')
    end
end
h=gscatter(double(catMethods),tmp.m,tmp.Size); legend off
for i=1:length(h)
    h(i).XData = h(i).XData + offset(i);
    h(i).Marker = S(i);
    h(i).MarkerSize= MSz;
    h(i).MarkerFaceColor = C{i};
    h(i).MarkerEdgeColor = 'k';%C(i) ;
    h(i).Color = 'k';%C(i);
end
t5.Title.String='Complexity = Power';
t5.Title.FontSize=26;
t5.YLim=[0 1.1];
t5.XLim=[0.5 14.5];
t5.XTick=1:14;
t5.XTickLabel=categories(catMethods);
t5.XTickLabelRotation=30;
t5.FontSize=26;
box on;

leg=legend(h,'Orientation', 'vertical');
leg.Layout.Tile = 'east';
leg.Interpreter = 'tex';
leg.Title.String="Size";
box on;
T1.XLabel.String="Method";T1.XLabel.FontWeight="bold";T1.XLabel.FontSize=26;
T1.YLabel.String="D_{method}"+newline;T1.YLabel.FontWeight="bold";T1.YLabel.FontSize=26;
