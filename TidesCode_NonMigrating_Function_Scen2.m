function [SeparationLT, PhaseError, AmplitudeError] = TidesCode_NonMigrating_Function_Scen2(Satellites)

% Solar Tides Computation - Migrating and Non-Migrating case.
% In this section we explore the solar tides in the SFC system by specifing 
% as set of tides and then plotting the individual modes. We sample
% the temperature function resulting from a sum of tide modes 
% and try to compute from the limited number of samples the amplitues and
% phases of the tide modes.
% In this version of the code we study the migrating and non-migrating modes.
% We included non-migrating tides of one earth rotation period and half-period.

% Define a set of parameters for the solar tides with respect to angular
% units of Local Time (LT) where there are 24 LT Hours over 2*pi radians

Mean = 1.0;               % Mean value which the tides oscillate about (K)
global r_lim
r_lim = 1.5*Mean;         % Upper limit of the r value in the polar plot    
Noise = 0.001;           % SNR = 30,20 and 10 dB or 10%, 1%, 0.1%, signal to "Kolmogorov Waves" ratio 
m_max = 100;              % m-number maximum of the "Kolmogorov waves"  
Mode = [1 2 3];           % List of migrating tide modes to simulate
SMode = [-3 -2 -1 0 1 2 3];    % List of non-Migrating tide modes  
Amplitude = [0.02 0.02 0.02; 0.02 0.02 0.02; 0.02 0.02 0.02; ... 
    0.02 0.02 0.02; 0.02 0.02 0.02; 0.02 0.02 0.02; 0.02 0.02 0.02]; 
                           % Amplitudes for each composite tide mode. 
                           % Each element Amplitude_ij means the amplitude
                           % of the i-non-migrating/j-migrating mode
Phase =[10 1 -1;  11 2 -3; 5 5 -2; ...
    -11 3 -2; -5 6 1; 6 3 2; -5 -4 2]; % Phase for each composite tide mode. 
                           % Each element Phase_ij means the phase of the
                           % i-non-migrating/j-migrating mode

LT = linspace(0, 24, 501); % Generate 501 linearly spaced points of Local 
                           % Times (LT) for plotting purposes, etc.
UT = linspace(0, 24, 501); % Generate 501 linearly spaced points of Universal 
                           % Times (UT) for plotting purposes, etc.
tide_names = cell(7,3);
tide_names(:,1) = {'DE3'; 'DE2'; 'DE1'; 'D0'; 'DW1'; 'DW2'; 'DW3'};
tide_names(:,2) = {'SE3'; 'SE2'; 'SE1'; 'S0'; 'SW1'; 'SW2'; 'SW3'};
tide_names(:,3) = {'TE3'; 'TE2'; 'TE1'; 'T0'; 'TW1'; 'TW2'; 'TW3'};


if ( (length(Mode)~=length(Amplitude(1,:)) ) || (length(Mode) ~= length(Phase(1,:))))
    err('The amplitude and phase matrix must be consistent with number of modes')
end

if ( (length(SMode) ~= length(Amplitude(:,1)) ) || (length(SMode) ~= length(Phase(:,1)) ) )
    err('The amplitude and phase matrix must be consistent with number of modes')
end

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This piece of code contains the fitting algorithm of the genearal problem
% of nom-migrating and migrating waves being sampled by satellites crossing
% the equator plane during a 24-hour period. The number of satellites,
% orbits and orbit period may be set up. We consider the satellites to be
% equally spaced in LT. The interval of time between two subsequent
% crossings over the equatorial plane is varied between t= 0 (synchronous case)
% to t = (12/Orbits)/Satellites (equally spaced assynchronous case). 
% Here the value t is equal for all satellites,  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = SolarTide(LT,UT, Mode,SMode, Amplitude,Phase,Mean); % Evaluate the temperature at the local times
KV_SMode = cat(2, [ (-m_max):(-4) ], [ 4:m_max ]);
KV_Mode = [1 2 3];
KV_Phase = zeros(length(KV_SMode),length(KV_Mode)); % Random phase of the Kolmogorov waves
for i = 1:length(KV_SMode)
    for k=1:length(KV_Mode)
        KV_Phase(i,k) = -12 + 24*rand();
    end
end
spT = KolmogorovWaves(LT,UT, Amplitude, m_max, Noise, KV_Phase); % high m-number waves (mod(m)>4) having Kolmogorov -5/3 Power Spectrum     


% Local time drift per day of satellite orbit 550 km altitude 70 deg inc
LTdrift = -2.54; % deg/day
Orbits = 15;     % orbits/day

SeparationLT = 0:0.001:12/Satellites; % The satellites are equally spaced in LT
SeparationUT = (12/Orbits)/Satellites; % Interval between equatorial crossings in UT

% Allocate memory for the conditioning matrix, and 
% error matrices
Acond = zeros(size(SeparationLT));
PhaseError = zeros(length(Mode)*length(SMode),length(SeparationLT));
AmplitudeError = zeros(size(PhaseError));

% Allocate memory for matrices representing the uncertainties estimated via
% covariance matrix
FitAmplitude = zeros(length(Mode)*length(SMode), length(SeparationLT));
FitPhase = zeros(length(Mode)*length(SMode), length(SeparationLT));
FitTemp = zeros(1, length(SeparationLT));

% Matrix  ( length(Mode).length(SMode) x 2 ) representing the ordinate pairs of the modes
ModeMat = zeros(length(Mode)*length(SMode),2);
for j=1:length(SMode)
    for i=1:length(Mode)
        ModeMat(i+length(Mode)*(j-1), 1) = SMode(j); 
        ModeMat(i+length(Mode)*(j-1), 2) = Mode(i)-SMode(j);
    end
end

% Reshape the amplitude and phase initial matrices, converting it into a
% 1-D array
AmpRef = reshape(Amplitude',[],1);
PhaseRef = reshape(Phase',[],1);

UTsat0 = zeros(2*Orbits,Satellites); % Allocate memory for initial starting UT matrix
for k=1:length(SeparationLT)
    LTsat0 = (0:(Satellites-1)) * SeparationLT(k); % Inital Starting LT of each satellite
    % Initial Starting UT time of each satellite
    for j = 1:2*Orbits
        UTsat0(j,:) = (j-1)*12/Orbits + (0:(Satellites-1)) * SeparationUT;
    end
    SUT = mod(reshape(UTsat0',1,[]),24); % put crossing times in 0 to 24 hour range and a 1D array.
    SLT = [LTsat0 LTsat0+12]' + (0:(Orbits-1)) .* (LTdrift/Orbits) * (24/360); % add drift to each orbit ascending and decending local times
    SLT = mod(reshape(SLT,1,[]),24); % put LT times in 0 to 24 hour range and a 1D array 
    Ts = zeros(1,length(SLT));
	spTs = zeros(1,length(SLT));
    for i=1:length(SLT)
        Ts(i) = SolarTide(SLT(i),SUT(i),Mode,SMode,Amplitude,Phase,Mean); % Temperature at each sample point
		spTs(i) = KolmogorovWaves(SLT(i),SUT(i), Amplitude, m_max, Noise, KV_Phase);
    end
    
    % Use the sampled values to find the amplitude and phase of tide modes

    % Create the sampling matrix in the order [Mean, Cos(), -Sin()]
    A = [ones(size(SLT)) ; cos(ModeMat(:,1) .* (pi/12*SLT) + ModeMat(:,2) .* (pi/12*SUT) ) ; sin(ModeMat(:,1) .* (pi/12*SLT) + ModeMat(:,2) .* (pi/12*SUT))]'; % MR: transpose - convert to (NDATA x NPARAM) matrix

    % find the condition of matrix A
    % The condition number associated with the linear equation Ax= b gives a 
    % bound on how inaccurate the solution x will be after approximation. 
    % one should think of the condition number as being (very roughly) the 
    % rate at which the solution, x, will change with respect to a change in b.
    % Thus, if the condition number is large, even a small error in b may cause 
    % a large error in x. On the other hand, if the condition number is small 
    % then the error in x will not be much bigger than the error in b. 
    Acond(k)= cond(A); % compute and save the conditon of the sampling matrix


    % we have the matrix equation samples = A * XY where samples and A are
    % potentially tall matrices or perhaps underdetermined. We use the
    % Moore-Penrose pseudoinverse of the A matrix to solve this equation.

    B = (A'*A)\(A'*(Ts'+spTs'));  % Compute weighted least squares.
    
    SMean = B(1);                 % pull out the mean from the solution
    X = B(2:1+max(size(ModeMat)));     % pull out X = Amplitude*Cos(Phase)
    Y = B(2+max(size(ModeMat)):end);   % pull out Y = Amplitude*Sin(Phase)

    SAmplitude = sqrt(X'.^2 + Y'.^2);      % Compute the sampled amplitude for each mode
    SPhase = 12/pi*(atan2(Y',X'));  % Compute the sample phase for each mode in units of LT
    
    % Save the input and output to see how close they are
    PhaseError(:,k) = 100*abs((PhaseRef-SPhase')./PhaseRef);
    AmplitudeError(:,k) = 100*abs((AmpRef-SAmplitude')./AmpRef);
    
    % Save the fitted parameters for plotting
    FitAmplitude(:,k) = SAmplitude';
    FitPhase(:,k) = SPhase';
    FitTemp(:,k) = SMean; 
end

make_ani_res = 0; % Set 1 in order to make this animation 
if (make_ani_res ~= 0) 
	filename_res = 'residuals_new_scen2.gif';
	h = figure('units','normalized','outerposition',[0 0 1 1]);
	s(1)= subplot(1,3,1,polaraxes);
	for i=1:1:length(UT)
		hold(s(1),'on')
		cla(s(1))
		polarplot(s(1),LT*pi/12,T(i,:),'LineWidth',2);
		polarplot(s(1),LT*pi/12,T(i,:)+spT(i,:),'.-','LineWidth',1);
		polarplot(s(1),LT*pi/12,ones(size(LT))*Mean,'--','LineWidth',1);
		hold(s(1),'off')
		s(1).ThetaTick = (0:23) * 360/24;        % set the polar plot ticks
		s(1).ThetaTickLabel = string(0:23);
		title(s(1),['Sampling of Solar Tides in Local Time Coordinates.' newline 'UT = ' num2str(UT(i))] );
		legend(s(1),'Solar Tide','Tide with noise','Mean Value');
		if (i == 1)
			s(3)= subplot(1,3,3); semilogy(SeparationLT,PhaseError,'linewidth',2);
			s(2)= subplot(1,3,2); semilogy(SeparationLT,AmplitudeError,'linewidth',2);
			title(s(2),'Amplitude Error Magnitude');
			title(s(3),'Phase Error Magnitude');
			xlabel(s(2),'Longitudinal Separation Between Satellites');
			xlabel(s(3),'Longitudinal Separation Between Satellites');
			ylabel(s(2),'Amplitude relative error (%)');
			ylabel(s(3),'Phase relative error (%)');
 
			legend(s(2),tide_names{1,1},tide_names{2,1}, tide_names{3,1} ...
				,tide_names{4,1}, tide_names{5,1}, tide_names{6,1} ... 
				,tide_names{7,1}, tide_names{1,2}, tide_names{2,2} ...
				,tide_names{3,2}, tide_names{4,2}, tide_names{5,2} ...
				,tide_names{6,2}, tide_names{7,2}, tide_names{1,3} ...
				,tide_names{2,3}, tide_names{3,3}, tide_names{4,3} ...
				,tide_names{5,3}, tide_names{6,3}, tide_names{7,3} ...
				,'Location','northeastoutside')
			%mytext = cell(size(Mode));
			%for k=1:length(Mode)
			%     mytext(k) = {['Mode-', num2str(Mode(k)), ' Mean=', num2str(Mean),' Ampl=', num2str(Amplitude(k)),' Phase=',num2str(Phase(k)),'LT']};
			%end
			%annotation('textbox',[0.6 0.52 0.38 0.06],'String',mytext,'FitBoxToText','on','BackgroundColor',[1 1 1]);
		end
		pause(0.1)
		% Capture the plot as an image 
		frame = getframe(h); 
		im = frame2im(frame); 
		[imind,cm] = rgb2ind(im,256); 
		% Write to the GIF File 
		if i == 1 
			imwrite(imind,cm,filename_res,'gif', 'Loopcount',inf); 
		else 
			imwrite(imind,cm,filename_res,'gif','DelayTime',0.1,'WriteMode','append'); 
		end
	end
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fuction definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a function to compute the Mean plus the Solar Tides (Migrating and Non-Migrating)
% Use matlab matrix operations to compute the sum of these cosine waves
% from the vector specifiers of Mode, Amplitude, Phase,and Mean
function T = SolarTide(LT,UT,Mode,SMode,Amplitude,Phase,Mean) 
    T = Mean + zeros(length(UT),length(LT));    
    for j = 1:length(UT)
        for l = 1:length(LT)
            for i = 1:length(SMode)
                for k=1:length(Mode)
                    T(j,l) = T(j,l) + Amplitude(i,k)* cos(2*pi/24*(SMode(i)*LT(l) + (Mode(k)-SMode(i))*UT(j) - Phase(i,k)));
                end
            end
        end
    end
    return     
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function to compute the "Kolmogorov Waves"
function spT = KolmogorovWaves(LT,UT, Amplitude, m_max, Noise, Phase) 
    a_signal = sqrt(sum(sum(Amplitude.^2)));
    m = 4:m_max;
    sum_m = sum(m.^(-5/3));
    a_0 = a_signal*sqrt(Noise/(6*sum_m));
    spT = zeros(length(UT),length(LT)); 
    SMode = cat(2, [ (-m_max):(-4) ], [ 4:m_max ]);
    Mode = [1 2 3];
    for j = 1:length(UT)
        for l = 1:length(LT)
            for i = 1:length(SMode)
                for k=1:length(Mode)             
                    spT(j,l) = spT(j,l) + abs(SMode(i))^(-5/6) * cos(2*pi/24*(SMode(i)*LT(l) + (Mode(k)-SMode(i))*UT(j) - Phase(i,k)));
                end
            end
        end
    end
    spT = a_0*spT;
    return     
end