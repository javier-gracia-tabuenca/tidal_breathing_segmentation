function [seg,volume,Vdrift] = cyclesAdvance(fs, volume, pars, varargin )
% [seg,volume,Vdrift] = cyclesAdvance(fs, volume, pars, varargin )
% Given a respiratory signal "volume" with sampling frequency "fs" detects 
% the begining of inspirations and expirations and divides the signal into
% respiratory cycles. On the process, it also detects and removes the
% volume base line. 
%
% This fuction is a secuence of 4 steps. First, estimate respiratory freq 
% based on the spectrogram. Second, simple estimation of begining of inspirations. 
% Third, -optional, recomended for IP signals- removes base line based on
% previous step's inspirations [1]. Forth, advance estimation of begining
% of inspiration and expirations [2]. 
%
% Further information below before each algorithm/step
%
% INPUTS:
% fs = sampling frequency
% volume =  respiratory volume signal
% pars = struct containing the parametes (explained below)
%
% OUTPUTS:
% seg.begIn = indices in "volume" signal where begining of inspirations occur
% seg.begEx = indices in "volume" signal where begining of expiration occur
%             It produces N begExp and N-1 endExp
% volume = volume signal with drift removed
% Vdrift = base line detected and removed
%
% VARARGIN
% 'plot' = plot final result
% 'plotDetails' = plot details on the 4 different steps 
% 'stopStep',n = stops the process after algorith number n
% 'timeThreshold' = 
%
% REFERENCES
%
% [1] Tidal breath analysis for infant pulmonary function testing
% J.H.T. Bates*,
% [2] Comparative Investigations of Algorithms for the Detection of
% Breaths in Newborns with Disturbed Respiratory Signals
% M. Schmidt, B
%
% Author: javier.gracia.tabuenca@gmail.com    Date: 12.12.2013
%



% default varargin
plotflag = '';
plotflag2 = '';
stopStep = 0;

for n = 1 : length(varargin)
	if strcmp(varargin{n}, 'plot')
		plotflag = 'plot';
	elseif strcmp(varargin{n}, 'plotDetails')
		plotflag2 = 'plot';
	elseif strcmp(varargin{n}, 'stopStep')
		n = n+1;
		stopStep = varargin{n};	
	elseif strcmp(varargin{n}, 'timeThreshold')
		n = n+1;
		timeThreshold = varargin{n};
	end
end


%default parameters
if isempty(pars)
	%pars step 1
	pars.estRF.Tresp_range_breathsMin = [6 50];
	%pars step 2
	pars.lowCutOff_Hz = 0.05;
	%pars step 3 (run if ~= 0)
	pars.baselineSpan_nPeaks = 5;
	%step 4
	pars.volThresholdIn_pMedian = 0.3;
	pars.volThresholdEx_pMedian = 0.3;
	pars.timeThresholdIn_pMedian = 0.3;
	pars.timeThresholdEx_pMedian = 0.3;
end

% Store the original volume signal
volume0 = volume;

%% STEP ONE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate the respiratory frequency from the spectrogram, as it is needed
% for the following step
%
% Parameters:
% pars.estRF.Tresp_range_breathsMin = [6 50];%find first resp rate between 6 and 50 breaths per minute
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
volume1 = volume0;

%default range [6 50] cylces per minute
[ Fresp1 ] = estimateRespFreq(fs, volume1, pars.estRF, plotflag2);
Tresp1 = 1/Fresp1;

%detailed plot
if strcmp(plotflag2,'plot')
	title(sprintf('STEP ONE: resp. freq. = %1.3f',Fresp1));
end


%% STEP TWO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate peaks/valleys from smoothed flow and volume;
% Removes volume drif using a low pass filter ; calcualtes flow ; smooths
% flow signal; find volume peaks as begExp ; finds valles between peak as
% begIns; 
%
% Quite rudimenatry method. Used only to detect the base line drift. 
% So that it can be removed in the following step as recomende in [1] 
%
% Parameters:
% pars.lowCutOff_Hz = 0.05; low pass filter cut off freq
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Remove low freq volume drift
n = 10; R = 60; Wst = pars.lowCutOff_Hz; % Order, stopband atten., stopband edge frequency (Hz)
[z,p,k] = cheby2(n,R,Wst/(fs/2), 'high');
[sos,g]=zp2sos(z,p,k);
h1=dfilt.df2sos(sos,g);
volume2 = filtfilthd(h1,volume1); % Zero-phase filtering

flow2 = diff(volume2)*fs;
flow2 = [flow2(1) ; flow2 ];

%smooth flow
span = Tresp1*fs*0.25;
flow2 = smooth(flow2,span);

%find begEx peaks over 0.2stds in the high passed volume
begEx2=find(diff(sign(flow2)) == -2);
begEx2(find(volume2(begEx2) < std(volume2)*0.2)) = [];

%find corresponding valleys
begIn2 = nan(length(begEx2)-1,1);
for i=1:length(begEx2)-1
	[~,mix] = min(volume1(begEx2(i):begEx2(i+1)));
	begIn2(i)= begEx2(i)+mix-1;
end

%re-find corresponding peaks
begEx2 = nan(length(begIn2)-1,1);
for i=1:length(begIn2)-1
	[~,mix] = max(volume1(begIn2(i):begIn2(i+1)));
	begEx2(i)= begIn2(i)+mix-1;
end


%detailed plot
if strcmp(plotflag2,'plot')
	flow0 = diff(volume0)*fs;
	flow0 = [flow0(1) ; flow0];
	t=getT(flow0,fs);
	figure;
	
	%volume
	h(1)=subplot(2,1,1);
	hold on
	plot(t,volume0,'b')
	plot(t,volume2,'c')
	plot(t(begEx2),volume0(begEx2),'.g')
	plot(t(begIn2),volume0(begIn2),'.r')
	grid on
	title('STEP TWO: Volume')
	xlabel('time [s]')
	ylabel('volume [a.u.]')
	legend({'original' 'smoothed' 'begExp' 'begIn'},'Location','eastoutside')
	%flow
	h(2)=subplot(2,1,2);
	hold on
	plot(t,flow0,'b')
	plot(t,flow2,'c')
	plot(t(begEx2),flow0(begEx2),'.g')
	plot(t(begIn2),flow0(begIn2),'.r')
	grid on
	title('Flow')
	xlabel('time [s]')
	ylabel('flow [a.u.]')	
	legend({'original' 'smoothed' 'begExp' 'begIn'},'Location','eastoutside')
	%
	linkaxes(h,'x');
	
end

% if enougth, leave here
if stopStep == 2	
	volume = cumtrapz(flow2);	
	Vdrift = [];
	seg.begIn = begIn2;
	seg.begEx = begEx2;
	return
end



volume3 = volume0;
begEx3 = begEx2;
begIn3 = begIn2;
Vdrift = [];
if pars.baselineSpan_nPeaks ~= 0
	%% STEP THREE : remove baseline on volume;
	% uses valles estimated in previous step to infer base line ; removes
	% base line; 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	n_point_mean = 2;

	x2 = [1:length(volume0)]';
	x1 = [1 ; begIn2 ; length(volume0) ];
	y1 = volume0(begIn2);
	y1 = [mean(y1(1:n_point_mean)) ; y1 ; mean(y1(end-n_point_mean:end)) ];
	%remove non unique
	[~,ix]=unique(x1);
	x1 = x1(ix);
	y1 = y1(ix);
	
	%Local regression using weighted linear least squares and a 1st degree polynomial model
	%A robust version of 'loess' that assigns lower weight to outliers in the regression.
	%The method assigns zero weight to data outside six mean absolute deviations.
	y1=smooth(y1,pars.baselineSpan_nPeaks,'rlowess')';
	
	%interpolate peaks
	Vdrift = interp1(x1,y1,x2,'spline');
	
	
	volume3 = volume0 - Vdrift;
	flow3 = diff(volume3)*fs;
	flow3 = [flow3(1) ; flow3 ];
	
	%didnt change
	begEx3 = begEx2;
	begIn3 = begIn2;
	
	%plot step 3
	if strcmp(plotflag2,'plot')
		t=getT(flow0,fs);
		figure;
		
		%volume
		h(1)=subplot(2,1,1);
		hold on
		plot(t,volume0,'c')
		plot(t,Vdrift,'r')
		plot(t,volume3,'b')
		grid on
		title('STEP TWO: Volume')
		xlabel('time [s]')
		ylabel('volume [a.u.]')
		legend({'no-baseline' 'original' 'baseline'},'Location','eastoutside')
		%flow
		h(2)=subplot(2,1,2);
		hold on
		plot(t,flow0,'k')
		plot(t,flow3,'b')
		grid on
		title('Flow')
		xlabel('time [s]')
		ylabel('flow [a.u.]')
		legend({'no-baseline' 'original'},'Location','eastoutside')
		%
		linkaxes(h,'x');
	end
	
	% if enougth, leave here
	if stopStep == 3
		volume = volume3;
		Vdrift = Vdrift;
		seg.begIn = begIn3;
		seg.begEx = begEx3;
		return
	end
	
end


%% STEP FOUR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Detects begining inspiration and expirations based on algorithm 4b [2].
% Detects all zero clossing in flow; rejects these zero-crossings which are
% under a time and volume threashold. Time and volume thresholds are
% respectively calcualted as 
% "pars.timeThresholdIn_pMedian/pars.timeThresholdIn_pMedian" 
% times the median of all inspiration/expiration periods and 
% "pars.volThresholdIn_pMedian/pars.volThresholdIn_pMedian" 
% times the median % of all inspiration/expiration volumes. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
volume4 = volume3;
flow4 = diff(volume4);
flow4 = [flow4(1) ; flow4];

% set to ins exp equal arrays
if begIn3(1)>begEx3(1) 	begEx3(1)=[];  end
l = min([length(begIn3) length(begEx3)]);


%Calculate volume aplitudes and durations for all expirations
tExs=zeros(1,l);
vExs=zeros(1,l);
for i=1:l-1
	tExs(i) = begIn3(i+1)-begEx3(i);
	vExs(i) = volume4(begIn3(i+1))-volume4(begEx3(i));
end
%Calculate volume aplitudes and durations for all inspirations
tIns=zeros(1,l);
vIns=zeros(1,l);
for i=1:l
	tIns(i) = begEx3(i)-begIn3(i);
	vIns(i) = volume4(begEx3(i))-volume4(begIn3(i));
	fIns(i) = max(flow4(begIn3(i):begEx3(i)));
end

% Thresholds
volThresholdIn = median(vIns)*pars.volThresholdIn_pMedian;
volThresholdEx = median(vExs)*pars.volThresholdEx_pMedian;

timeThresholdIn= median(tIns)*pars.timeThresholdIn_pMedian;
timeThresholdEx = median(tExs)*pars.timeThresholdEx_pMedian;

%detect zeroCross
[zcross,zcrossD]=zeroCrossAdvance(flow4);
zcrossN=zcross(find(zcrossD==-1));

%check wrong zeroCrossings 
begIn4=[];
begEx4=[];
for i=2:length(zcross)-1
	
	if zcrossD(i)>0 %temporal begIn
		%criteria
		if volThresholdIn < volume4(zcross(i+1))-volume4(zcross(i)) && ...
				timeThresholdIn < zcross(i+1)-zcross(i)
			
			begIn4=[begIn4 zcross(i)];
		end
	else
		%criteria
		if (volThresholdEx > volume4(zcross(i-1))-volume4(zcross(i)) && ...
				timeThresholdEx < zcross(i)-zcross(i-1))
			
			begEx4=[begEx4 zcross(i)];
		end
	end
end


%remove duplicates
rix = find(begEx4<begIn4(1));
begEx4(rix)=[];

rix = find(begEx4>begIn4(end));
begEx4(rix)=[];

i=1;
while i+1 <= length(begIn4) && i <= length(begEx4)
	[~,ix]=find( begIn4 > begIn4(i) & begIn4 < begEx4(i));
	begIn4(ix)=[];
	
	[~,ix]=find( begEx4 > begEx4(i) & begEx4 < begIn4(i+1));
	begEx4(ix)=[];
	
	i=i+1;
end
% begIn5 one more length than begEx5
begIn4=begIn4(1:i);
begEx4=begEx4(1:i-1);

% !Correction :
% algorithm 4b [2]
% Thresholds in the oposite direction left to rigth for the epiratory

begInb=begIn4;
begExb=begEx4;

begIn4=begIn4(1);
begEx4=[];

for i=1:length(begExb)
	
	if (volThresholdEx > volume4(begInb(i+1))-volume4(begExb(i)) && ...
			timeThresholdEx < begInb(i+1)-begExb(i))
		
		begEx4=[begEx4 begExb(i)];
		begIn4=[begIn4 begInb(i+1)];
	end
	
end
	
%plot step 3
if strcmp(plotflag2,'plot')
	t=getT(flow4,fs);
	figure
	
	h(1)=subplot(2,1,1);
	hold on
	plot(t,volume4,'b')
	plot(t(zcross),volume4(zcross),'xr')
	plot(t(zcrossN),volume4(zcrossN),'xg')
	plot(t(begIn4),volume4(begIn4),'or')
	plot(t(begEx4),volume4(begEx4),'og')
	grid on
	title('STEP FOUR: Volume')
	xlabel('time [s]')
	ylabel('volume [a.u.]')
	legend({'volume' '+0cross' '-0cross' 'begIn' 'begExp' },'Location','eastoutside')
	
	h(2)=subplot(2,1,2);
	hold on
	plot(t,flow4,'b')
	plot(t(zcross),flow4(zcross),'xr')
	plot(t(zcrossN),flow4(zcrossN),'xg')
	plot(t(begIn4),flow4(begIn4),'or')
	plot(t(begEx4),flow4(begEx4),'og')
	plot(t,flow4,'b')
	grid on
	title('STEP FOUR: Volume')
	xlabel('time [s]')
	ylabel('flow [a.u.]')
	legend({'flow' '+0cross' '-0cross' 'begIn' 'begExp' },'Location','eastoutside')
	
	linkaxes(h,'x')
	
end
	

begEx=begEx4;
begIn=begIn4;
flow=flow4;
volume=volume4;

%% Create return struct
seg.begIn = begIn;
seg.begEx = begEx;

%% plot
if strcmp(plotflag,'plot') || strcmp(plotflag2,'plot')
	t=getT(flow,fs);
	figure;
	
	%volume
	h(1)=subplot(2,1,1);
	hold on
	plot(t,volume,'b')
	plot(t(begEx),volume(begEx),'.r')
	plot(t(begIn),volume(begIn),'.m')
	grid on
	title('Volume')
	xlabel('time [s]')
	ylabel('volume [a.u.]')
	legend({'volume' 'begIn' 'begExp' },'Location','eastoutside')
	%flow
	h(2)=subplot(2,1,2);
	hold on
	plot(t,flow,'b')
	plot(t(begEx),flow(begEx),'.r')
	plot(t(begIn),flow(begIn),'.m')
	grid on
	title('Flow')
	xlabel('time [s]')
	ylabel('flow [a.u.]')
	legend({'flow' 'begIn' 'begExp' },'Location','eastoutside')
	%
	linkaxes(h,'x');
end





end


%zero cross detector
function [zcross,zcrossD]=zeroCrossAdvance(x)

% buld array zcross with zero crossinindexes,
zcross = find(abs(diff(sign(x))) ==  2);
% and array with crsossign direction negatives when crosing from +to-
zcrossN = find(diff(sign(x)) == -2);
[~,locsN] = ismember(zcrossN,zcross);
zcrossD = ones(size(zcross));
zcrossD(locsN)=-1;
%correction for low sampling freq: between two points of crossing, get the closest to 0
for i=1:length(zcross)
	[~,ixm]=min(abs([ x(zcross(i)) x(zcross(i)+1) ]));
	if ixm==2
		zcross(i)=zcross(i)+1;
	end
end

end






