clear all
close all
clc

[fPath, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(fPath)
addpath([fPath,'\Functions'])
addpath([fPath,'\Model'])
%% Get the raw audio signal from the synchronized recordings

dataDir = {
    [fPath,'\Data\R50B_0+1_20221119$134500'];
    [fPath,'\Data\R32A_0+1_20221119$134500'];
    [fPath,'\Data\R50A_0+1_20221119$134500'];
    [fPath,'\Data\R32C_0+1_20221119$134500'];
    };

data_Test = cell(height(dataDir),1);

for k = 1:height(dataDir)
    dataset = fullfile(dataDir{k,1}); 
    ads = audioDatastore(dataset,'IncludeSubfolders',true);
    
    for i = 1:height(ads.Files)
        [audioIn,dsInfo] = read(ads);
        Fs = dsInfo.SampleRate;
         
    end

    data_Test{k,1} = audioIn(:,1); %Process channel 1
end


%% Process the recordings using the detection module

% Call the detection module
[Detect_Loc, label_Pred] = Detection_Module(data_Test, Fs);

%% Match the detected calls using the matching module

% Parameter for matching
Temp=6;                 % Field Temperature
v = 331.3+0.606*Temp;   % Speed of sound (meters/second)
t_avg = 0.15;           % Average Duration of a call (seconds)
t_min = 0.05;           % Minimum Duration of a call (seconds)
% c = 5;                % Maximum Duration of a call, defined as c*t_avg (integer)
D_max = 113;            % Maximum distance between any pairs of microphones (meters)
GPS_Error = 3;          % GPS error for microphones' locations measurement (meters)
epsilon = 0.1;          % Min duration between possible windows (seconds)
t_max = t_avg*5;
D_max = D_max + 2*GPS_Error;    % Calibrate GPS Error
min_seperation_length = D_max/v*0.7;  % Minimum duration between separable sets of calls

%-------------------------------%

% Call the matching module
[Matched_Pred_Calls,Matched_Pred_Vec] = Matching_Module(Detect_Loc,data_Test,Fs,Temp,t_avg,t_min,t_max,D_max,min_seperation_length,epsilon);

%% Visualize one of the set of clips containing matched calls 

close all

choice = 1;
for k = 1:height(Matched_Pred_Calls) 
    close all
    % If choice ~= 1, sequentially display all sets
    if choice ~= 1
        break
    end
    a = k;
    lt = Matched_Pred_Calls{a,"Start Time"}-1;
    ht = Matched_Pred_Calls{a,"End Time"}+1;
    
    f = figure('Name',['Visualize set of clips #',num2str(a),'/',num2str(height(Matched_Pred_Calls)),': ',num2str(lt),':',num2str(ht)],'WindowState','maximized');
    pause(0.1)
    
    for i = 1:height(data_Test)
        subplot(height(data_Test),1,i,'Parent',f) 
        
        pred = label_Pred{i,1}(1,round(lt*Fs):round(ht*Fs));    
        intersect_P = double(Matched_Pred_Vec(round(lt*Fs):round(ht*Fs),1))*3+1;
        hold on
        spectrogram(data_Test{i,1}(round(lt*Fs):round(ht*Fs),1),Fs*0.04,round(0.8*Fs*0.04),1000:1:3500,Fs,'yaxis')
        [s,~,t,p] = spectrogram(data_Test{i,1}(round(lt*Fs):round(ht*Fs),1),Fs*0.04,round(0.8*Fs*0.04),1000:1:3500,Fs,'yaxis');
        zm = max(max(abs(s)));
        
        z = zeros(1,width(t));
        z(1,:) = zm;
        f_p = pred(1,round(t*Fs));
        f_i_p = intersect_P(round(t*Fs),1);
    
        index = f_p == 2;
        f_p(index) = 3;
        f_p(~index) = 1;
        
        hold on
        plot3(t,f_p,z,'Color',[1, 0, 0, 0.7],'LineWidth',0.01)
        plot3(t,f_i_p,z,'Color',[0, 0, 0, 0.7],'LineWidth',0.01)
    
        if ht-lt > 1
            xticks(0:(ht-lt)/10:(ht-lt));
            xticklabels(num2str((lt:(ht-lt)/10:ht)'));
        else
            xticks((0:(ht-lt)/10:(ht-lt))*1000);
            xticklabels(num2str((lt:(ht-lt)/10:ht)'));
            xlabel('Time (s)')
        end
        title(['Audio Sample ',num2str(lt),':',num2str(ht),', Rec#',num2str(i)])
    
        legend('Predicted-Label','Matched-Pred','Location','northeastoutside')
        hold off
    
    end
    % choice = inputdlg([ 'Displaying set ',num2str(a),'/',num2str(height(Matched_Pred_Calls)),newline,...
    %     'Enter 1 to continue displaying sets, any other values to skip']);
    % choice = str2double(choice{1,1});
    if a < height(Matched_Pred_Calls)   % There are remaining sets not displayed
        choice = questdlg(['Displaying set ',num2str(a),'/',num2str(height(Matched_Pred_Calls)),newline,...
                           'Press continue to continue displaying sets, press skip to skip'], ...
                           'Choice','continue','skip','continue');
        if strcmp(choice,'continue')
            choice = 1;
        else
            choice = 0;
        end
    else
        f = msgbox(['Displaying set ',num2str(a),'/',num2str(height(Matched_Pred_Calls))]);
        choice = 0;
    end
end
%% Find the TDOA between the matched calls in the sets of matched clips
%  using the TDOA module

% Call the TDOA module
TimeLag_Table = TDOA_Module(data_Test,Fs,D_max,v,Detect_Loc,Matched_Pred_Calls);

%% Localize the source using TDOA with localization module

resolution = 1;

range_offset = 150;

Mic = [
352253.61	3734845.58  % WNE
% 352171.16	3734846.85  % WNW
352250.54	3734766.87  % WSE
% 352169.98	3734770.34  % WSW

% 352432.42	3734843.92  % MNE
352353.68	3734845.14  % MNW
% 352431.19	3734764.08  % MSE
352352.45	3734765.3   % MSW

% 352429.65	3734664.28  % SNE
% 352349.02	3734663.31  % SNW
% 352430.27	3734584.41  % SSE
% 352349.67	3734585.66  % SSW
];

% Call the localization module
[cum_result_matrix,range,localization_log, Pred_Loc] = Localization_Module(TimeLag_Table.Lag_Table{:,:}, TimeLag_Table.Corr_Table{:,:}, Mic, Temp, resolution, range_offset);

mic = [352253.61	3734845.58  % WNE
352171.16	3734846.85  % WNW
352250.54	3734766.87  % WSE
352169.98	3734770.34  % WSW

352432.42	3734843.92  % MNE
352353.68	3734845.14  % MNW
352431.19	3734764.08  % MSE
352352.45	3734765.3   % MSW

352429.65	3734664.28  % SNE
352349.02	3734663.31  % SNW
352430.27	3734584.41  % SSE
352349.67	3734585.66  % SSW
];

figure('Name','Map');
hold on

plot(mic(:,1),mic(:,2),'s',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0.5,0.5,0.5], ...
    'DisplayName','Mic Loc')

%----------------- Plot true loc and pred loc-------------------%
True_Loc = [                   
    352339.63	3734793.45
];
plot(True_Loc(1),True_Loc(2),'o',...
    'LineWidth',2,...
    'MarkerSize',5,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',[0.5,0.5,0.5], ...
    'DisplayName','True Loc with 10 m circle')
plot(Pred_Loc(:,1),Pred_Loc(:,2),'o',...
    'LineWidth',2,...
    'MarkerSize',5,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5], ...
    'DisplayName','Pred Loc')
viscircles(True_Loc,10,'Color','r','LineWidth',1)
%-------------------------------------------------------------%

axis([range(1,1) range(2,1) range(1,2) range(2,2)]);
axis equal
grid on
legend
xlabel('East (m)')
ylabel('North (m)')

