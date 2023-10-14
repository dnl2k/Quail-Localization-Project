
function TimeLag_Table = TDOA_Module(data_Test,Fs,D_max,v,Detect_Loc,Matched_Pred_Calls)
    [fPath, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
    addpath(fPath)
    TimeLag_Table = [];
    TimeLag_Corr_Table = [];
    for a = 1:height(Matched_Pred_Calls)
        offset = 0.2; %Seconds
        R_threshold = 0.1; % Threshold as percent difference between 2 highest peaks of cross correlation coefficient 
        lt = Matched_Pred_Calls{a,"Start Time"}-offset;
        ht = Matched_Pred_Calls{a,"End Time"}+offset;
        
        % Find First Detection
        First_Group_ST = [];
        for j = 1:height(data_Test)
            if numel(Matched_Pred_Calls{a,j}{1}) >0
                First_Group_ST = [First_Group_ST; j, Detect_Loc{j,1}{Matched_Pred_Calls{a,j}{1}(1),4}];
            end
        end
        First_Detection_Index = find(First_Group_ST(:,2) == min(First_Group_ST(:,2)));
        
        dummy = First_Group_ST(:,1)';
        dummy(dummy == First_Group_ST(First_Detection_Index(1),1)) = [];
        dummy = [First_Group_ST(First_Detection_Index(1),1), dummy]; %First Detection signal will be reference microphone value
        
        % Get the raw audio signals
        refsig = data_Test{dummy(1,1),1}(round(Fs*lt):round(Fs*ht),1);
        
        sig = zeros(height(refsig),width(dummy)-1);
        for i = 2:width(dummy)
            sig(:,i-1) = data_Test{dummy(1,i),1}(round(Fs*lt):round(Fs*ht),1);
        end
        
        % Bandpass the signals
        refsig = bandpass(refsig,[1000 3500],Fs); 
        for i = 1:width(sig)
            sig(:,i) = bandpass(sig(:,i),[1000 3500],Fs);
        end
        
        % Remove correlated noise using canned Wiener Filter function
        [refsig,~] = WienerNoiseReduction(refsig,Fs,Fs*offset);
        for i = 1:width(sig)
            [sig(:,i),~] = WienerNoiseReduction(sig(:,i),Fs,Fs*offset);
        end
        
        % Find the TDOA using normalized cross-correlation
        tau_est = zeros(1,width(sig));        
        tau_corr = zeros(1,width(sig));
        for i = 1:width(sig)
            [r,lags] = xcorr(refsig,sig(:,i),'normalized',round(D_max/v*Fs));
            index = find(abs(r) == max(abs(r(lags < 0 & lags > -round(D_max/v*Fs)))));
            tau_est(:,i) = abs(lags(index)/Fs);
            tau_corr(:,i) = r(index);
        end
    
        % Return the TDOA and corresponding coefficient
        tau = NaN(1,12);
        tau(1,dummy(1,1)) = 0;
        tau_coefficient = NaN(1,12);
        tau_coefficient(1,dummy(1,1)) = 0;
        
        for i = 2:width(dummy)
            if tau_est(1,i-1) > 0
                tau(1,dummy(1,i)) = tau_est(1,i-1);
                tau_coefficient(1,dummy(1,i)) = tau_corr(1,i-1);
            end
        end
                
        TimeLag = array2table(tau,'VariableNames',["Mic 1","Mic 2","Mic 3","Mic 4",...
            "Mic 5","Mic 6","Mic 7","Mic 8",...
            "Mic 9","Mic 10","Mic 11","Mic 12"]);
        TimeLag_Corr = array2table(tau_coefficient,'VariableNames',["M1 Corr","M2 Corr","M3 Corr","M4 Corr",...
            "M5 Corr","M6 Corr","M7 Corr","M8 Corr",...
            "M9 Corr","M10 Corr","M11 Corr","M12 Corr"]);

        TimeLag_Table = [TimeLag_Table;TimeLag(1,1:height(data_Test))];
        TimeLag_Corr_Table = [TimeLag_Corr_Table;TimeLag_Corr(1,1:height(data_Test))];
     
    end
    dummy = TimeLag_Table;
    TimeLag_Table = [];
    TimeLag_Table.Lag_Table = dummy;
    TimeLag_Table.Corr_Table = TimeLag_Corr_Table;
end