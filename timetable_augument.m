function [X_lagged,Y_lagged,XY_lagged]=timetable_augument(X,Y,Delta)
%----------------------------------------------------------------------maca
%------------------------------------------------------------------<header>
% Function to augment timetables X and Y with lagged variables.
% Inputs:
%   X - Original timetable for independent variables
%   Y - Original timetable for dependent variables
%   Delta - Structure containing parameters for lagging
%-----------------------------------------------------------------<\header>

% Initialize variables and adjust lag parameters based on Delta
[obs,var]=size(X);
dLower=Delta.Lower; % Lower bound for lags
dUpper=Delta.Upper; % Upper bound for lags
if min(dLower)<0
    lagYFlag=1; % Flag indicating need to lag Y
    lagYAmount=-min(Delta.Lower);% Amount of lag for Y
    % Adjust upper and lower bounds for laggingdUpper=dUpper-min(Delta.Lower);
    dLower=dLower-min(Delta.Lower);
    dUpper=dUpper-min(Delta.Lower);
else
    lagYFlag=0;
    lagYAmount=0;
end
max_lag=max(Delta.max_lag,max(dUpper));% Maximum lag to be applied
if obs> max_lag
    % Initialize lagged timetables for X and Y
    X_lagged=timetable('Size',[obs,var*max_lag],'VariableTypes',repmat("doublenan",[1,var*max_lag]),'RowTimes',X.Properties.RowTimes);
    Y_lagged=timetable('Size',[obs,size(Y,2)],'VariableTypes',repmat("doublenan",[1,size(Y,2)]),'RowTimes',Y.Properties.RowTimes);
    Y_lagged.Properties.VariableNames=Y.Properties.VariableNames;
else
    % Warning if not enough observations to apply lag
    warning(("The minimum number of observations is max_lag" + newline...
        + "Not enough observations: Aborting."));
    return
end
% Preparing array to store variable names for lagged timetables
names_array=cell(var*max_lag,1);
Duration=mode(diff(X.Properties.RowTimes));% Common duration between time steps
intrvl=[0;find(diff(X.Properties.RowTimes)~=Duration);obs];% Identify intervals of consistent duration
k=1;
for l=1:length(intrvl)-1
    if (intrvl(l+1)-intrvl(l))>=max_lag
        % Extract parts of X and Y for current interval
        X_Temp=X(intrvl(l)+1:intrvl(l+1),:);
        Y_Temp=Y(intrvl(l)+1:intrvl(l+1),:);
        % Apply lag to Y if needed
        if lagYFlag
            Y_partial=lag(Y_Temp,lagYAmount);
        else
            Y_partial=Y_Temp;
        end
        Y_lagged(intrvl(l)+1:intrvl(l+1),:)=Y_partial;
        % Apply lags to X and fill in X_lagged
        for i=1:max_lag
            X_partial=lag(X_Temp,i-1);
            for j=1:var
                if i>dLower(j)&&i<=dUpper(j)
                    idx=max_lag*(j-1)+i;
                    X_lagged(intrvl(l)+1:intrvl(l+1),idx)=X_partial(:,j);
                    idx_names=max_lag*(j-1)+i;
                    if k==1
                        if lagYAmount==0
                            if i>1
                                names_array{idx_names}=X.Properties.VariableNames{j} + " (t-" + num2str(i) + ")";
                                
                            else
                                names_array{idx_names}=X.Properties.VariableNames{j} + " (t)";
                            end
                        else
                            if i>lagYAmount+1
                                names_array{idx_names}=X.Properties.VariableNames{j} + " (t-" + num2str(i-lagYAmount-1) + ")";
                                
                            elseif i==lagYAmount+1
                                names_array{idx_names}=X.Properties.VariableNames{j} + " (t)";
                            else
                                names_array{idx_names}=X.Properties.VariableNames{j} + " (t+" + num2str(lagYAmount-i+1) + ")";
                            end
                        end
                        
                    end
                else
                    idx_names=max_lag*(j-1)+i;
                    if k==1
                        if lagYAmount==0
                            if i>1
                                names_array{idx_names}=X.Properties.VariableNames{j} + " (t-" + num2str(i-1) + ")";
                                
                            else
                                names_array{idx_names}=X.Properties.VariableNames{j} + " (t)";
                            end
                        else
                            if i>lagYAmount+1
                                names_array{idx_names}=X.Properties.VariableNames{j} + " (t-" + num2str(i-lagYAmount-1) + ")";
                                
                            elseif i==lagYAmount+1
                                names_array{idx_names}=X.Properties.VariableNames{j} + " (t)";
                            else
                                names_array{idx_names}=X.Properties.VariableNames{j} + " (t+" + num2str(lagYAmount-i+1) + ")";
                            end
                        end
                        
                    end
                end
            end
            
            
        end
        
        k=0;% Flag to indicate completion of one interval
    end
    
end

X_lagged.Properties.VariableNames=cellstr(names_array); % Update variable names in X_lagged
X_lagged=X_lagged(max_lag+1:end,:); % Trim initial rows with NAN data because of the lagging
Y_lagged=Y_lagged(max_lag+1:end,:); % Same trimming for Y_lagged
% Apply filters to remove rows with missing or infinite values 
% #MIGHT NEED AN ERROR IF INF SHOW UP - Why would it show up??
filter0=~(any(ismissing(X_lagged),1));
filter1=~(any(ismissing(X_lagged(:,filter0)),2)|any(ismissing(Y_lagged),2));
filter2=~(any(isinf(X_lagged{:,:}),2)|any(isinf(Y_lagged{:,:}),2));
filter3=filter1&filter2;
X_lagged=X_lagged(filter3,filter0);
Y_lagged=Y_lagged(filter3,:);
% Combine lagged X and Y into one timetable
XY_lagged=[X_lagged,Y_lagged];
end
