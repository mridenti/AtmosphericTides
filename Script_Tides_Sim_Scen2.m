close all;
clear all;

Satellites = 4; %varied from 4 to 6

SIM_ITER = 500; %500, default value

% First Iteration
[SeparationLT, PhaseError, AmplitudeError] = TidesCode_NonMigrating_Function_Scen2(Satellites);
PhaseErrorBuffer = PhaseError;
AmplitudeErrorBuffer = AmplitudeError;

style= {'-b','-k','-r', '-g', '-y', '-c', '-m',...
        '--b','--k','--r', '--g', '--y', '--c', '--m', ...
        ':b',':k',':r', ':g', ':y', ':c', ':m'};

tide_names(:,1) = {'DE3'; 'DE2'; 'DE1'; 'D0'; 'DW1'; 'DW2'; 'DW3'};
tide_names(:,2) = {'SE3'; 'SE2'; 'SE1'; 'S0'; 'SW1'; 'SW2'; 'SW3'};
tide_names(:,3) = {'TE3'; 'TE2'; 'TE1'; 'T0'; 'TW1'; 'TW2'; 'TW3'};

for i=2:SIM_ITER
    [SeparationLT, PhaseError, AmplitudeError] = TidesCode_NonMigrating_Function_Scen2(Satellites);
    PhaseErrorBuffer = PhaseErrorBuffer + PhaseError;
    AmplitudeErrorBuffer = AmplitudeErrorBuffer + AmplitudeError;
end

PhaseError = PhaseErrorBuffer/SIM_ITER;
AmplitudeError = AmplitudeErrorBuffer/SIM_ITER;
dim_matrix = size(AmplitudeError);

s(1) = subplot(1,2,1); 
hold on
for k=1:dim_matrix(1)
    semilogy(SeparationLT,AmplitudeError(k,:),style{k},'linewidth',2);    
end
s(2) = subplot(1,2,2); 
hold on
for k=1:dim_matrix(1)
    semilogy(SeparationLT,PhaseError(k,:),style{k},'linewidth',2);    
end

xlabel(s(1),'Longitudinal Separation Between Satellites');
xlabel(s(2),'Longitudinal Separation Between Satellites');
ylabel(s(1),'Amplitude relative error (%)');
ylabel(s(2),'Phase relative error (%)');

legend(s(2),tide_names{1,1},tide_names{2,1}, tide_names{3,1} ...
            ,tide_names{4,1}, tide_names{5,1}, tide_names{6,1} ... 
            ,tide_names{7,1}, tide_names{1,2}, tide_names{2,2} ...
            ,tide_names{3,2}, tide_names{4,2}, tide_names{5,2} ...
            ,tide_names{6,2}, tide_names{7,2}, tide_names{1,3} ...
            ,tide_names{2,3}, tide_names{3,3}, tide_names{4,3} ...
            ,tide_names{5,3}, tide_names{6,3}, tide_names{7,3} ...
            ,'Location','northeastoutside')
        