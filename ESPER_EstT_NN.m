function [Estimates]= ...
    ESPER_EstT_NN(DesiredVariables,OutputCoordinates,PredictorMeasurements,PredictorTypes, ... % Required inputs
            varargin)                                     % Optional inputs
% This function is not extensively documented because it is only intended
% to be called from within TRACEv1.  It will estimate temperature from
% salinity and coordinate information.
    % Input Parameter Key: 
    % 1. Salinity

VerboseTF=false;

% *************************************************************************
% Parsing inputs, setting defaults, and sanity checking inputs.
%
% Starting timer
tic

% Verifying required inputs are provided
if nargin<4; error('ESPER_EstT called with too few input arguments.'); end

% Checking whether specific equations are specified.
a=strcmpi(varargin,'Equations');
if any(a)
    Equations=varargin{1,logical([0 a(1:end-1)])};
else
    Equations=1;
end
% Making [] argument for Equations equivalent to no argument.
if isempty(Equations);Equations=1:4; end
% Making 0 argument for Equations equivalent to no argument.
if Equations==0; Equations=1:4; end

% Checking for PerKgSwTF input and setting default if not given
a=strcmpi(varargin,'PerKgSwTF');
if any(a)
    PerKgSwTF=varargin{1,logical([0 a(1:end-1)])};
else
    PerKgSwTF=true;
end

% ESPER_EstT requires non-NaN coordinates to provide an estimate.  This step
% eliminates NaN coordinate combinations prior to estimation.  NaN
% estimates will be returned for these coordinates.
NaNGridCoords=max(isnan(OutputCoordinates),[],2);

% Doing a size check for the coordinates.
if ~(size(OutputCoordinates,2)==3)
    error('OutputCoordinates has either too many or two few columns.  This version only allows 3 columns with the first being longitude (deg E), the second being latitude (deg N), and the third being depth (m).');
end

% Figuring out how many estimates are required
n=sum(~NaNGridCoords);
e=size(Equations,2);
p=size(DesiredVariables,2);

% Checking for common missing data indicator flags and warning if any are
% found.  Consider adding your commonly used flags here.
if max(ismember([-999 -9 -1*10^20],PredictorMeasurements))==1
    warning('ESPER_EstT: A common non-NaN missing data indicator (e.g. -999, -9, -1e20) was detected in the input measurements provided.  Missing data should be replaced with NaNs.  Otherwise, ESPER_EstT will interpret your inputs at face value and give terrible estimates.')
end

% Book-keeping with coordinate inputs and adjusting negative longitudes.
C=OutputCoordinates(~NaNGridCoords,:);
C(C(:,1)>360,1)=rem(C(C(:,1)>360,1),360);
C(C(:,1)<0,1)=rem(C(C(:,1)<0,1,1),360)+360;

% Throwing an error if latitudes are out of the expected range.
if max(abs(C(:,2)))>90
    error('ESPER_EstT: A latitude >90 degrees (N or S) has been detected.  Verify latitude is in the 2nd colum of the coordinate input.');
end

% % Preparing full PredictorMeasurement uncertainty grid
% DefaultUncertainties=diag([1 1 0.02 0.02 0.02 0.01]);
% DefaultUAll=zeros(size(PredictorMeasurements,1),6);
% DefaultUAll(:,PredictorTypes)=PredictorMeasurements* ...
%     DefaultUncertainties(PredictorTypes,PredictorTypes);                   % Setting multiplicative default uncertainties for P, N, O2, and Si.
% DefaultUAll(:,ismember(PredictorTypes,[1 2]))=0.003;                       % Then setting additive default uncertainties for T and S.
% DefaultUAll=DefaultUAll(~NaNGridCoords,:);
% if UseDefaultUncertainties==false
%     InputUAll=zeros(size(PredictorMeasurements));
%     InputUAll(:,PredictorTypes)=InputU;
%     InputUAll=max(cat(3,InputUAll, DefaultUAll),[],3);                     % Overriding user provided uncertainties that are smaller than the (minimum) default uncertainties
% else
%     InputUAll=DefaultUAll;
% end  

% Making sure all provided predictors are identified.
if ~(size(PredictorTypes,2)==size(PredictorMeasurements,2))
    error('The PredictorTypes input does not have the same number of columns as the Measurements input.  This means it is unclear which measurement is in which column.');
end

% Putting all provided measurement inputs in standard order
MAll=NaN(n,6);                                                             % PredictorMeasurements
% UAll=NaN(n,6);                                                             % Uncertainties
MAll(:,PredictorTypes)=PredictorMeasurements(~NaNGridCoords,:);            % Reordering and limiting to viable coordinates
% UAll(:,PredictorTypes)=InputUAll(:,PredictorTypes);                        % This was already limited to viable coordinates for later use.
YouHaveBeenWarnedCanth=false;                                              % Calculating Canth is slow, so this flag is used to make sure it only happens once.

% Beginning the iterations through the requested properties
for PIter=1:p                                                              
    Property=DesiredVariables(1,PIter);
    % Specifying which variables are required for each property estimate.
    % Each row is for a different property.
    NeededForProperty=[1 2 3 6 5
        1 2 3 6 5
        1 2 3 6 5
        1 2 3 6 5
        1 2 3 6 5
        1 2 3 6 5
        1 2 3 6 5
        1 2 3 6 5];
    % Which of the 5 predictors are required for the equations specified?
    % Each row is for a different equation.
    VarVec=logical([1 1 0 0 0]); % depth sal and temp in NN
    NeedVars=any(VarVec(Equations,:),1);
    HaveVars=false(1,6);  HaveVars(1,PredictorTypes+1)=1; % Adding one because depth is provided
    HaveVars(1,1)=1; % This is depth in this version
    
    % Temperature is required if O2 is used as a predictor (to convert to
    % AOU)... equation 15 is therefore only valid under some circumstances.
    if ismember(6,PredictorTypes)==1; NeedVars(1,2)=true; end
    
    % Ignoring variables that aren't used for the current property
    HaveVarsSub=HaveVars(1,NeededForProperty(Property,:));

    % Temperature is required if measurements are provided in molar units
    % or if CO2sys will be used to adjust pH values
    if PerKgSwTF==false; NeedVars(1,2)=1; end 

    % Making sure all needed variables are present
    if ~all(HaveVarsSub(1,NeedVars)) && VerboseTF==true  % We are assuming only power users turn VerboseTF off, and hopefully they understand the function call
        disp('Warning: One or more regression equations for the current property require one or more input parameters that are either not provided or not labeled correctly. These equations will return NaN for all estimates.  All 16 equations are used by default unless specific equations are specified.  Temperature is also required when density or carbonate system calculations are called.'); % They should already be NaN
    end
    
    % Limiting measurements to the subset needed for this property estimate
    M=MAll(:,NeededForProperty(Property,:));
%     U=UAll(:,NeededForProperty(Property,:));
%     DefaultU=DefaultUAll(:,NeededForProperty(Property,:));

    % Checking to see whether temperature is needed, and subbing in
    % potential temp if it is.  As of the version 3 update the  routines
    % now request temperature, but convert it to potential temperature
    % before using it for predictions.
    if NeedVars(1,2)
        M(:,2)=sw_ptmp(M(:,1),M(:,2),sw_pres(C(:,3),C(:,2)),0);
    end
    % Checking to see whether O2 is needed. Defining AOU and subbing in for
    % O2 if yes (see above).
    if any(ismember(NeededForProperty(Property,:),6))
        M(:,4)=sw_satO2(M(:,1),M(:,2))*44.64-M(:,4);
    end
    % Converting units to molality if they are provided as molarity.
    if PerKgSwTF==false
        densities=sw_dens(M(:,1),M(:,2),sw_pres(C(:,3),C(:,2)))/1000;
        M(:,3)=M(:,3)./densities;
        M(:,4)=M(:,4)./densities;
        M(:,5)=M(:,5)./densities;
    end

    % *********************************************************************
    % Beginning treatment of inputs and calculations
    if     Property==1; VName='Temperature';
    else; error('A property identifier >1 or <1 was supplied, but this routine only has 1 possible property estimate.  The property identifier is the first input.')
    end
    
    % Loading the data, with an error message if not found
    FN='Polys.mat';
    % Making sure you downloaded the needed file and put it somewhere it
    % can be found
    if exist(FN,'file')<2; error('ESPER_EstT could not find the file(s) needed to run.  These mandatory file(s) should be distributed from the same website as ESPER_EstT.  Contact the corresponding author if you cannot find it there.  If you do have it then make sure all of the contents of the ESPER.zip extract are on the MATLAB path or in the active directory.  This will require adding several subfolders for ESPER.'); end
    L=load(FN);
  
    
    %Preallocating for speed
    OutputEstimates=NaN(size(OutputCoordinates,1),e);                      % Using size instead of n since we want to preallocate for NaN coordinate locations as well   
    Est=NaN(n,e);

    % Some of the equations use depth [m] as a predictor so this appends it
    % as a predictor.  This, so far, is unique to TRACE.
    M=horzcat(C(:,3),M);  


    % Disambiguation:
    % Eq... a counter for which equation ESPER_EstT is on
    % e... the number of equations that will be used
    % Equation... the specific equation number currently being used
    % Equations... the user provided list of equations that will be used.
    for Eq=1:e                                                             % Iterating over the (up to 16) equations
        Equation=Equations(1,Eq);
        NeedVarsNow=any(VarVec(Equation,:),1);
        if all(HaveVarsSub(1,NeedVarsNow))                                 % Skipping if we don't have what we need
            P=horzcat(cosd(C(:,1)-20),sind(C(:,1)-20),C(:,2),M(:,VarVec(Equations(Eq),:)));
            for Net=1:4 % A committee of 4 neural networks is used.
                % Separate neural networks are used for the Arctic/Atlantic and
                % the rest of the ocean (more info below).  This functional
                % form of ESPER_NN is designed to avoid calling the nerual
                % network toolbox.  It will require the (many) ESPER_NN
                % functions to be on the MATLAB path.
                FunctionOther=str2func(horzcat('TRACE_EstT_',VName,'_',num2str(Equation),'_Other_',num2str(Net))); % Basically an eval statement.  Sorry Matlab power users.
                FunctionAtl=str2func(horzcat('TRACE_EstT_',VName,'_',num2str(Equation),'_Atl_',num2str(Net)));
                EstAtl(:,Eq,Net)=FunctionAtl(P')';
                EstOther(:,Eq,Net)=FunctionOther(P')';
            end
        end
    end
    % Averaging across neural network committee members
    EstAtl=nanmean(EstAtl,3);
    EstOther=nanmean(EstOther,3);

    % We do not want information to propagate across the Panama Canal (for
    % instance), so data is carved into two segments... the Atlantic/Arctic
    % (AA) and everything else.
    AAInds=or(inpolygon(C(:,1),C(:,2),L.Polys.LNAPoly(:,1),L.Polys.LNAPoly(:,2)), ...
        or(inpolygon(C(:,1),C(:,2),L.Polys.LSAPoly(:,1),L.Polys.LSAPoly(:,2)), ...
        or(inpolygon(C(:,1),C(:,2),L.Polys.LNAPolyExtra(:,1),L.Polys.LNAPolyExtra(:,2)), ...
        or(inpolygon(C(:,1),C(:,2),L.Polys.LSAPolyExtra(:,1),L.Polys.LSAPolyExtra(:,2)), ...
        inpolygon(C(:,1),C(:,2),L.Polys.LNOPoly(:,1),L.Polys.LNOPoly(:,2))))));
    % We'd like a smooth transition in the Bering Strait and in the South
    % Atlantic.  This linearly interpolates between the networks at their
    % boundaries.
    BeringInds=inpolygon(C(:,1),C(:,2),L.Polys.Bering(:,1),L.Polys.Bering(:,2));
    SoAtlInds=(C(:,1)>290 | C(:,1)<20) & (C(:,2)>-44) & (C(:,2)<-34);
    SoAfrInds=(C(:,1)>19 & C(:,1)<27) & (C(:,2)>-44) & (C(:,2)<-34);
    Est=EstOther;                                                          % Pulling out Indo-Pacific estimates
    Est(AAInds,:)=EstAtl(AAInds,:,:);                                      % Pulling out all other estimates
    Est(BeringInds,:)=EstAtl(BeringInds,:).* ...                           % Blending estimates in the Bering Strait vicinity
        repmat((C(BeringInds,2)-62.5)/(7.5),[1,size(Equations,2)]) ...
        +EstOther(BeringInds,:).* ...
        repmat((70-C(BeringInds,2))/(7.5),[1,size(Equations,2)]);
    Est(SoAtlInds,:)=EstAtl(SoAtlInds,:).* ...                             % Blending estimates in the Southern Atlantic to Southern Ocean transition
        repmat((C(SoAtlInds,2)+44)/(10),[1,size(Equations,2)])+ ...
        EstOther(SoAtlInds,:).* ...
        repmat((-34-C(SoAtlInds,2))/(10),[1,size(Equations,2)]);
    Est(SoAfrInds,:)=Est(SoAfrInds,:).* ...                                % Blending estimates South of South Africa toward Indo-Pac Estimates
        repmat((27-C(SoAfrInds,1))/(8),[1,size(Equations,2)])+ ...
        EstOther(SoAfrInds,:).* ...
        repmat((C(SoAfrInds,1)-19)/(8),[1,size(Equations,2)]);
    
    % Filling the outputs with the estimates at viable locations.
    OutputEstimates(~NaNGridCoords,:)=Est;
    Estimates.(VName)=OutputEstimates;
end