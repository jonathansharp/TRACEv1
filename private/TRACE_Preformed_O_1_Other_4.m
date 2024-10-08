function [Y,Xf,Af] = TRACE_Preformed_O_1_Other_4(X,~,~)
%TRACE_PREFORMED_O_1_OTHER_4 neural network simulation function.
%
% Auto-generated by MATLAB, 21-Aug-2024 10:13:38.
% 
% [Y] = TRACE_Preformed_O_1_Other_4(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 6xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 1xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-0.999961923064171;-0.999961923064171;-77.5;0;30.1649227142334;-2.9811577796936];
x1_step1.gain = [1.00003807838574;1.00003807838574;0.0119760479041916;0.000363636363636364;0.22152510192902;0.0532708976942295];
x1_step1.ymin = -1;

% Layer 1
b1 = [1.4875548313409043821;0.77439383114511706196;1.004204222089130516;1.0152010168769602938;-16.368418628166452322];
IW1_1 = [-0.055244142402583316931 0.012558021706086197264 0.097850298277788474999 -0.011374988678288651747 0.4680798171897025961 1.0349251435358477469;-0.23488634600808294683 0.08420877853339148178 -0.48387282518472707338 0.10838510573508906842 0.3111516029144504869 0.94945877143425805311;0.3930502646120718202 -0.17515499893683378341 0.70077299625484756884 -0.28432979603562319459 -2.3388763759341619775 1.2095459724284383096;-0.16797433734462877131 0.062113778435089889163 0.51543468858710683733 0.01150986087347444578 0.11566910312873940336 -0.34282388860367724037;-0.037184517205050816868 0.017534613809449393906 0.17847965578469096593 -16.21607541045525025 0.011702090690938608472 -0.20453298409532780222];

% Layer 2
b2 = [3.0176806280586525588;-0.78004440622203308209;2.1227875524182127087;2.7247684045270759157;3.49804901471492391;-2.2232974958213826966;0.56722898884941419517;-2.0682394881427574873;1.3951842755521557393;-1.9019390075863396028;-2.2738992756976386111;-1.9474467244751338946;4.757829702745664413;6.1768549738175435593;2.6214058655140561882];
LW2_1 = [-2.2969331032722517527 0.97424108396413089395 -0.12220006285262743906 -1.2824620300775411152 -0.077659633708075354042;-5.6899602433165723525 4.9925256630954537584 3.5876452507343210563 7.2133101464635549505 2.6260520714920474461;1.1044754938478071882 0.28438437748788436554 0.21692677667738369252 -2.2278911974039652399 0.85311756091580992489;-2.3925186568732588199 0.56121241056339654651 -0.28264627251813506659 -1.0766200991570797907 -0.211601014285572403;2.4070005634670121886 -0.42263695103079912396 -0.25710773281705845417 -4.7578273042201226772 0.66871978191012737547;-1.3290475673912280907 -0.62641277455265620411 0.021147603598116179041 2.6238794130327915965 -1.9894332745800482609;-1.2211477759588373893 1.7228124351149263216 1.3067662043275238837 1.5747655607608372907 1.4280789936669429085;9.4093436405413299894 -1.1430188882170078113 -1.3620097810701592422 -5.3874919171921824557 -1.0926359020721374549;0.076600018927166044214 1.0915084202492972842 0.83280016673835954855 -0.37486947523899899481 1.2730990192083395662;10.278795567652904808 -1.4216353562127082011 -1.633153425499981104 -6.0176237094242068792 -0.94837427619485914487;11.108470316385334442 -4.4755394913563968373 -2.7822475625366509533 -3.9320810054915944498 -0.62913473124210728127;12.844205777758201492 -5.2964034615641173431 -3.3986738357222083806 -5.2153060925018124294 -0.59528851334943577722;4.0094434223800572425 0.60473916052787468622 3.7000653674334502696 -8.918188506580202457 -6.0162677436108742413;-11.859073821088129463 3.7565530114587812882 -3.6238835216511753323 -0.73716111422842678369 -19.065199130439850705;-9.0884357447671817454 3.4692453945881505106 2.1018077944749151698 2.4242679127082675805 0.55413256439856661828];

% Layer 3
b3 = -0.92530024283140532848;
LW3_2 = [-3.8991514145501424693 -0.65125679679662451527 14.266844540328534663 4.5266123465499887857 -3.4059680765447915007 0.99783293758406654828 8.536831787053669629 1.9030580066191584709 -14.726975766549690761 -1.7769195809022086952 -9.0894298415916310319 4.1710533702434551273 -0.09457250962725285981 -0.16112660801754316586 -5.3403736607517737767];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00612908122568115;
y1_step1.xoffset = 89.6956371703064;

% ===== SIMULATION ========

% Format Input Arguments
isCellX = iscell(X);
if ~isCellX
  X = {X};
end

% Dimensions
TS = size(X,2); % timesteps
if ~isempty(X)
  Q = size(X{1},2); % samples/series
else
  Q = 0;
end

% Allocate Outputs
Y = cell(1,TS);

% Time loop
for ts=1:TS

    % Input 1
    Xp1 = mapminmax_apply(X{1,ts},x1_step1);
    
    % Layer 1
    a1 = tansig_apply(repmat(b1,1,Q) + IW1_1*Xp1);
    
    % Layer 2
    a2 = tansig_apply(repmat(b2,1,Q) + LW2_1*a1);
    
    % Layer 3
    a3 = repmat(b3,1,Q) + LW3_2*a2;
    
    % Output 1
    Y{1,ts} = mapminmax_reverse(a3,y1_step1);
end

% Final Delay States
Xf = cell(1,0);
Af = cell(3,0);

% Format Output Arguments
if ~isCellX
  Y = cell2mat(Y);
end
end

% ===== MODULE FUNCTIONS ========

% Map Minimum and Maximum Input Processing Function
function y = mapminmax_apply(x,settings)
  y = bsxfun(@minus,x,settings.xoffset);
  y = bsxfun(@times,y,settings.gain);
  y = bsxfun(@plus,y,settings.ymin);
end

% Sigmoid Symmetric Transfer Function
function a = tansig_apply(n,~)
  a = 2 ./ (1 + exp(-2*n)) - 1;
end

% Map Minimum and Maximum Output Reverse-Processing Function
function x = mapminmax_reverse(y,settings)
  x = bsxfun(@minus,y,settings.ymin);
  x = bsxfun(@rdivide,x,settings.gain);
  x = bsxfun(@plus,x,settings.xoffset);
end
