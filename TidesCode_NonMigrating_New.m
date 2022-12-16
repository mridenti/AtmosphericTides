% Solar Tides Computation - Migrating and Non-Migrating case.
% In this section we explore the solar tides in the SFC system by specifing 
% as set of tides and then plotting the individual modes. We sample
% the temperature function resulting from a sum of tide modes 
% and try to compute from the limited number of samples the amplitues and
% phases of the tide modes.
% In this version of the code we study the migrating and non-migrating modes.
% We included non-migrating tides of one earth rotation period, half-period and quarter-period.

% Define a set of parameters for the solar tides with respect to angular
% units of Local Time (LT) where there are 24 LT Hours over 2*pi radians

Mean = 1.0;               % Mean value which the tides oscillate about (K)
global r_lim
r_lim = 1.5*Mean;           % Upper limit of the r value in the polar plot  
Noise = 0.1;           % SNR = 30,20 and 10 dB or 10%, 1%, 0.1%, signal to "Kolmogorov Waves" ratio
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
% This first animation shows each mode separately evolving in UT time. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

make_ani_A = 0; % Set 1 in order to make this animation 
if (make_ani_A ~= 0)
scale = Mean/(max(max(Amplitude))*2) ;
for k = 1:length(Mode) % Loop over each tide modes and put in next window
    filename_modes = ['modes_mig_nonmig_new_n' num2str(k) '.gif'];
    h = figure('units','normalized','outerposition',[0 0 1 1]);
    for i=1:1:length(UT)
        for j=1:length(SMode)
           subplot(ceil(sqrt(length(SMode))),ceil(sqrt(length(SMode))),j);
           LTplot(LT,SolarTide(LT, UT(i), Mode(k), SMode(j), Amplitude(j,k)*scale, Phase(j,k), Mean),'LineWidth',1);
           hold on
           LTplot(LT,ones(size(LT))*Mean,'--','LineWidth',1) 
           hold off
           title({strcat('Tide Mode-','(',num2str(Mode(k)),',',num2str(SMode(j)),')', ...
             ' - ', tide_names{j,k} ,' UT= ', num2str(UT(i))); strcat('Mean=',num2str(Mean), ...
               ' Ampl=',num2str(Amplitude(j,k)),' Phase=',num2str(Phase(j,k)),'LT')},'FontSize',6);
        end    
        % Capture the plot as an image 
        frame = getframe(h); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 
        % Write to the GIF File 
        if i == 1 
            imwrite(imind,cm,filename_modes,'gif', 'Loopcount',inf); 
        else 
            imwrite(imind,cm,filename_modes,'gif','DelayTime',0.1,'WriteMode','append'); 
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This second animation shows the sum of the modes evolving in UT time.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

make_ani_B = 0; % Set 1 in order to make this animation 
if (make_ani_B ~= 0)
close all
filename_sim = 'simulated_tides_new.gif';
h = figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:2:length(UT)
    LTplot(LT,SolarTide(LT, UT(i), Mode, SMode, Amplitude, Phase, Mean),'LineWidth',1);
    hold on
    LTplot(LT,ones(size(LT))*Mean,'--','LineWidth',1) 
    hold off
    title({strcat('UT= ', num2str(UT(i))); strcat('Mean=',num2str(Mean))},'FontSize',6);
 
  % Capture the plot as an image 
  frame = getframe(h); 
  im = frame2im(frame); 
  [imind,cm] = rgb2ind(im,256); 
  % Write to the GIF File 
  if i == 1 
       imwrite(imind,cm,filename_sim,'gif', 'Loopcount',inf); 
  else 
       imwrite(imind,cm,filename_sim,'gif','DelayTime',0.1,'WriteMode','append'); 
  end
end
end

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
spT = KolmogorovWaves(LT,UT, Amplitude, m_max, Noise); % high m-number waves (mod(m)>4) having Kolmogorov -5/3 Power Spectrum   

% Local time drift per day of satellite orbit 550 km altitude 70 deg inc
LTdrift = -2.54; % deg/day
Orbits = 15;     % orbits/day
Satellites = 4;  % Number of satellite in constelation

SeparationLT = 12/Satellites; % The satellites are equally spaced in LT
SeparationUT = 0:0.05:(12/Orbits)/Satellites; % Interval between equatorial crossings in UT

% Allocate memory the  conditioning matrix, and 
% error matrices
Acond = zeros(size(SeparationUT));
PhaseError = zeros(length(Mode)*length(SMode),length(SeparationUT));
AmplitudeError = zeros(size(PhaseError));

% Allocate memory for matrices representing the uncertainties estimated via
% covariance matrix
FitAmplitude = zeros(length(Mode)*length(SMode), length(SeparationUT));
FitPhase = zeros(length(Mode)*length(SMode), length(SeparationUT));
FitTemp = zeros(1, length(SeparationUT));

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
for k=1:length(SeparationUT)
    LTsat0 = (0:(Satellites-1)) * SeparationLT; % Inital Starting LT of each satellite
    % Initial Starting UT time of each satellite
    for j = 1:2*Orbits
        UTsat0(j,:) = (j-1)*12/Orbits + (0:(Satellites-1)) * SeparationUT(k);
    end
    SUT = mod(reshape(UTsat0',1,[]),24); % put crossing times in 0 to 24 hour range and a 1D array.
    SLT = [LTsat0 LTsat0+12]' + (0:(Orbits-1)) .* (LTdrift/Orbits) * (24/360); % add drift to each orbit ascending and decending local times
    SLT = mod(reshape(SLT,1,[]),24); % put LT times in 0 to 24 hour range and a 1D array 
    Ts = zeros(1,length(SLT));
    spTs = zeros(1,length(SLT));
    for i=1:length(SLT)
        Ts(i) = SolarTide(SLT(i),SUT(i),Mode,SMode,Amplitude,Phase,Mean); % Temperature at each sample point
        spTs(i) = spT(k,i);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This animation illustrates roughly the time when the crossings occur for 
% each UT separation. The temperature profile is also shown and one can see
% its variation over time. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
make_ani_C = 0; % Set 1 in order to make this animation 
if (make_ani_C ~= 0)
close all
filename_cro = 'simulated_crossings_new.gif';
h = figure('units','normalized','outerposition',[0 0 1 1]); % for plotting each sampling scenerio
for i=1:1:length(UT)
    for k=1:length(SeparationUT)
        % Plot the sum of all the modes and the sample points
        % create a bunch of subplots based on the loop count
        % Initial Starting UT of each satellite
        % Inital Starting LT of each satellite
        LTsat0 = (0:(Satellites-1)) * SeparationLT;
        % Initial Starting UT of each satellite
        for j = 1:2*Orbits
             UTsat0(j,:) = (j-1)*12/Orbits + (0:(Satellites-1)) * SeparationUT(k);
        end
        SUT = mod(reshape(UTsat0',1,[]),24); % put crossing times in 0 to 24 hour range and a 1D array
        SLT = [LTsat0 LTsat0+12]' + (0:(Orbits-1)) .* (LTdrift/Orbits) * (24/360); % add drift to each orbit ascending and decending local times
        SLT = mod(reshape(SLT,1,[]),24); % put results in 0 to 24 hour range and a 1D array. 
        %SLT = [ 0 12 1 13 2 14 3 15]; % Local Times to Sample the temperature at (SLT)
        Ts = zeros(1,length(SLT));
        for j=1:length(SLT)
              Ts(j) = SolarTide(SLT(j),SUT(j),Mode,SMode,Amplitude,Phase,Mean); % Temperature at each sample point
        end
        Tother = normrnd(0,Noise,size(T));  
        Tsother = normrnd(0,Noise,size(Ts));% T_other noise like terms at sample times
        rows = floor(sqrt(length(SeparationUT)));
        colls = ceil(sqrt(length(SeparationUT)));
        subplot(rows,colls,k);
        LTplot(LT,T(i,:),'LineWidth',2);
        hold on
        LTplot(LT,T(i,:)+Tother(i,:),'.-','LineWidth',1);
        LTplot(LT,ones(size(LT))*Mean,'--','LineWidth',1) 
        for j = 1:length(SLT)
            if (SUT(j) <= UT(i)+UT(2) && SUT(j) >= UT(i)-UT(2))
                LTplot(SLT(j),Ts(j)+Tsother(j),'*','LineWidth',2) 
            end
        end
        hold off
        % if k==1; legend('Solar Tide','Tide with noise','Mean Value','Sample Points'); end; % only one legend
        title(['Separation ' num2str(SeparationUT(k)) 'UT' newline 'UT: ' num2str(UT(i))])
    end
    pause(0.1)
    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if i == 1 
        imwrite(imind,cm,filename_cro,'gif', 'Loopcount',inf); 
    else 
        imwrite(imind,cm,filename_cro,'gif','DelayTime',0.1,'WriteMode','append'); 
    end
end
end

close all
% Create output graph of condition number 
figure
plot(SeparationUT,Acond, '-o', 'linewidth',2)
title({'Sampling Matrix Condition Number'; ['Constelation of ' num2str(Satellites) ' satellites']},'FontSize',18);
ylabel('Condition','FontSize',18);
xlabel('Local Time Separation Between Satellites','FontSize',18);
pause()

make_ani_res = 0; % Set 1 in order to make this animation 
if (make_ani_res ~= 0) 
filename_res = 'residuals_new.gif';
h = figure('units','normalized','outerposition',[0 0 1 1]);
s(1)= subplot(1,3,1,polaraxes);
for i=1:1:length(UT)
    hold(s(1),'on')
    cla(s(1))
    polarplot(s(1),LT*pi/12,T(i,:),'LineWidth',2);
    polarplot(s(1),LT*pi/12,T(i,:)+Tother(i,:),'.-','LineWidth',1);
    polarplot(s(1),LT*pi/12,ones(size(LT))*Mean,'--','LineWidth',1);
    hold(s(1),'off')
    s(1).ThetaTick = (0:23) * 360/24;        % set the polar plot ticks
    s(1).ThetaTickLabel = string(0:23);
    title(s(1),['Sampling of Solar Tides in Local Time Coordinates.' newline 'UT = ' num2str(UT(i))] );
    legend(s(1),'Solar Tide','Tide with noise','Mean Value');
    if (i == 1)
        s(3)= subplot(1,3,3); semilogy(SeparationUT,PhaseError,'linewidth',2);
        s(2)= subplot(1,3,2); semilogy(SeparationUT,AmplitudeError,'linewidth',2);
        title(s(2),'Amplitude Error Magnitude');
        title(s(3),'Phase Error Magnitude');
        xlabel(s(2),'Universal Time Separation Between Satellites');
        xlabel(s(3),'Universal Time Separation Between Satellites');
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

% Now see how the fitting compares with the simulated curve
k = 1; % k = 1 tests the syncronous case
ReshAmplitude = reshape(FitAmplitude(:,k)', [length(Mode), length(SMode)])';
ReshPhase = reshape(FitPhase(:,k), [length(Mode), length(SMode)])';

make_ani_fit = 1; % Set 1 in order to make this animation 
if (make_ani_fit ~= 0)
close all
filename_fit = 'fitted_tides_four_satellites_4dB.gif';
h = figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:2:length(UT)
    sim_T = SolarTide(LT, UT(i), Mode, SMode, Amplitude, Phase, Mean) ...
        + KolmogorovWaves(LT,UT(i), Amplitude, m_max, Noise);
    LTplot(LT,sim_T,'k','LineWidth',1);
    hold on
    LTplot(LT,ones(size(LT))*Mean,'--','LineWidth',1) 
    LTplot(LT,SolarTide(LT, UT(i), Mode, SMode, ReshAmplitude, ReshPhase, FitTemp(k)),'--b','LineWidth',2);
    hold off
    legend('Simulated','Mean Temperature','Fitted curve');
    title({strcat('UT= ', num2str(UT(i))); strcat('Mean=',num2str(Mean))},'FontSize',6);
 
  % Capture the plot as an image 
  frame = getframe(h); 
  im = frame2im(frame); 
  [imind,cm] = rgb2ind(im,256); 
  % Write to the GIF File 
  if i == 1 
       imwrite(imind,cm,filename_fit,'gif', 'Loopcount',inf); 
  else 
       imwrite(imind,cm,filename_fit,'gif','DelayTime',0.1,'WriteMode','append'); 
  end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fuction definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a function to compute the Mean plus the Solar Tides (Migrating and Non-Migrating)
% Use matlab matrix operations to compute the sum of these cosine waves
% from the vector specifiers of Mode, Amplitude, Phase,and Mean
function T = SolarTide(LT,UT,Mode,SMode,Amplitude, Phase,Mean) 
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

% Create a function to compute the "Kolmogorov Waves"
function spT = KolmogorovWaves(LT,UT, Amplitude, m_max, Noise) 
    a_signal = sqrt(sum(sum(Amplitude.^2)));
    m = 4:m_max;
    sum_m = sum(m.^(-5/3));
    a_0 = a_signal*sqrt(Noise/(6*sum_m));
    spT = zeros(length(UT),length(LT));    
    SMode = cat(2, [ (-m_max):(-4) ], [ 4:m_max ]);
    Mode = [1 2 3];
    Phase = zeros(length(SMode),length(Mode));
    for i = 1:length(SMode)
        for k=1:length(Mode)
            Phase(i,k) = -12 + 24*rand();
        end
    end
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

% Auxiliary functions
function LTplot(LT, r, varargin)
%
% The inputs are
%    LT: array of local time values (0 to 24 hr) 
%    r:  radius vector
%    varargin: additional arguments to add onto the polarplot function 

 global r_lim

 polarplot(LT*2*pi/24,r,varargin{:});
 rlim([0 r_lim]);
 ax = gca;
 ax.ThetaTick = (0:23) * 360/24;        % set the polar plot ticks
 ax.ThetaTickLabel = string(0:23);      % set the polar plot labels
end
