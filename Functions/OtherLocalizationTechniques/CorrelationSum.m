function [Pred_Loc] = CorrelationSum(data_Test,Fs, Detect_Loc,Matched_Pred_Calls, Mic, Temp)
    % This function performs Correlation Sum localization methods, written 
    % based on algorithm from: 
    % Acoustic localization of antbirds in a Mexican rainforest using a wireless sensor network. 
    % Collier et al. 2010.
    %
    % Input: 
    % 1. data_Test: the recordings from each recorder
    % 2. Detect_Loc: table containing details about detected calls for each
    % recorder.
    % 3. Matched_Pred_Calls: table with details about the window containing
    % matched calls
    % 4. Mic: recorders' locations
    % 5. Temp: air temperature
    %
    % Output:
    % 1. Pred_Loc: list of the coordinates of the localized calls.

    Pred_Loc = [];
    for a =1:height(Matched_Pred_Calls) 
        v=331.3+0.606*Temp; % Speed of sound in dry air
        
        % Initialize the cross correlation coefficients
        C = cell(4,4); 
        C_hat = cell(4,4);
        C_lags = cell(4,4);
        
        data_Test_bp = data_Test;
        % Bandpass the tracks
        for i = 1:height(data_Test_bp)
            data_Test_bp{i,1} = bandpass(data_Test_bp{i,1},[1000 3500],Fs);
        end
        
        % Find First Detection as first key signal
        First_Group_ST = [];
        for j = 1:height(data_Test)
            if numel(Matched_Pred_Calls{a,j}{1}) >0
                First_Group_ST = [First_Group_ST; j, Detect_Loc{j,1}{Matched_Pred_Calls{a,j}{1}(1),4},Detect_Loc{j,1}{Matched_Pred_Calls{a,j}{1}(1),5}];
            end
        end
        First_Detection_Index = find(First_Group_ST(:,2) == min(First_Group_ST(:,2)));
        
        dummy = First_Group_ST(:,1)';
        dummy(dummy == First_Group_ST(First_Detection_Index(1),1)) = [];
        dummy = [First_Group_ST(First_Detection_Index(1),1), dummy]; %First Detection signal will be reference microphone value
        
        % Get signal
        lt = First_Group_ST(First_Group_ST(:,1) == dummy(1,1),2);
        ht = First_Group_ST(First_Group_ST(:,1) == dummy(1,1),3);
        key_sig = data_Test_bp{dummy(1,1),1}(round(Fs*lt):round(Fs*ht),1);
        
        
        k = dummy(1,1);
        t_k = lt;       % Start time of key signal
        s_k = ht-lt;    % Duration of the key signal
        
        T = dummy;      % All tracks with detected calls
        
        while width(T) > 1
            T(T== k) = [];      % Step 1: Remove k from T
            
            C_max = -99;
            l_max = -99;
            i_max = k;
            for i = 1:width(T)  % Step 2: Loop through every track in T
                % Calculate the parameter for X_i (recordings) segments
                d_ki = sqrt(sum((Mic(k,:)-Mic(T(i),:)).^2));
                lt = t_k - d_ki/v;
                ht = lt + s_k + 2*(d_ki/v);
                
                
                % Extract X_i
                X_i = data_Test_bp{T(i),1}(round(Fs*lt):round(Fs*ht),1);

                % Zero pad the key signal
                dummy = zeros(height(X_i),1);
                dummy(round(Fs*(d_ki/v))+(1:height(key_sig)),1) = key_sig; 
                
                % Compute the cross-correlation
                [C{T(i),k},C_lags{T(i),k}] = xcorr(dummy,X_i);
                C_hat{T(i),k} = abs(hilbert(C{T(i),k}));
            
                [M,I] = max(C{T(i),k});
                if M > C_max
                    C_max = M;
                    l_max = C_lags{T(i),k}(I)/Fs;
                    i_max = T(i);
                end
            end
            k = i_max;
            t_k = t_k + l_max;
        end
        
        % Map the cross correlation coefficient to the lattice
        
        range_offset = 100;
        resolution = 1;
        range=[min(Mic)-range_offset;max(Mic)+range_offset];
        [X, Y] = meshgrid(range(1,1):resolution:range(2,1), range(1,2):resolution:range(2,2));
        S = zeros(size(X));
        
        for i = 1:4
            for j = 1:4
                close
                if isempty(C{i,j})
                    continue
                end
                % Map the lattice's lag values
                dummy = (sqrt((X-Mic(i,1)).^2+(Y-Mic(i,2)).^2)-sqrt((X-Mic(j,1)).^2+(Y-Mic(j,2)).^2))/v;    
                C_lags_min = min(C_lags{i,j})/Fs;
                C_lags_max = max(C_lags{i,j})/Fs;
                dummy(dummy > C_lags_max | dummy < C_lags_min) = nan;
                imagesc(dummy);
                set(gca,'YDir','normal')
                colorbar
                
                % With the lattice's lag values, map the corresponding cross-corr
                % coefficient
                dummy = round(dummy*Fs);                        
                dummy_unique = unique(dummy);
                dict = containers.Map(C_lags{i,j},C_hat{i,j});
                for k = 1:height(dummy_unique)
                    dummy(dummy == dummy_unique(k)) = dict(-dummy_unique(k));
                end
        
                % Accumulate S
                S = S + dummy;
            end
        end
        
        % Map the cross-correlation coefficient to the higher resolution lattice 
        % centered around the previous step localization results
        
        index = find(S == max(S,[],'all'));
        
        Pred_Loc = [Pred_Loc;[(X(index)) (Y(index))]];
    end
end

