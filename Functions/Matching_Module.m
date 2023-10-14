
 
function [Matched_Pred_Calls,Matched_Pred_Vec] = Matching_Module(Detect_Loc,data_Test,Fs,Temp,t_avg,t_min,t_max,D_max,min_seperation_length,epsilon)
    %-------- Find the window to match the Predicted calls
    Matched_Pred_Calls_Times = find_Suitable_Calls(Detect_Loc,Temp,t_avg,t_min,t_max,D_max);
    
    %-------- Generate the Predicted mask using the window information
    Mask_Pred = generate_Label_Mask(Matched_Pred_Calls_Times,data_Test,Fs);
    
    %--------- Combine the Predicted window across the locations
    logical_Mask_Pred = generate_Logical_Label_Mask(Mask_Pred);
    [Intersect_Pred, Intersect_Pred_Vec] = generateIntersectCall(logical_Mask_Pred,Fs);
    
    %-------- Matched the Pred calls using intersected Window
    Matched_Pred_Calls = generate_Matched_Window(Intersect_Pred,Detect_Loc,data_Test);
    Matched_Pred_Calls = trim_Detection(Matched_Pred_Calls,data_Test);    % Trim detection without 3+ detections across 4 locations
    
    %------- Fix >2 calls group to a single window
    Matched_Pred_Calls = seperate_Calls(Matched_Pred_Calls,Detect_Loc,t_avg,min_seperation_length,epsilon,data_Test,D_max,Temp);
    
    %-------- Generate mask for visualization
    Matched_Pred_Vec = generate_Mask_Visualization(Matched_Pred_Calls,Intersect_Pred_Vec,Fs);
end

% Functions

function [Detected_Label,logical_label_Test_I] = generateIntersectCall(logical_label_Test,Fs)
    % This function take the logical arrays and sampling frequency and
    % return the regions with possible matched calls in the format of
    % logical vector and table

    logical_label_Test_I = double(logical_label_Test{1,1}) + double(logical_label_Test{2,1});
    for i = 3:height(logical_label_Test)
        logical_label_Test_I = logical_label_Test_I + double(logical_label_Test{i,1});    
    end
    index = logical_label_Test_I >= 3;
    logical_label_Test_I(index,:) = 1;
    logical_label_Test_I(~index,:) = 0;
    
    dummy = double(logical_label_Test_I);
    index = dummy(:,1) == 1;
    dummy(index,:) = 2;
    dummy(~index,:) = 1;
    
    rising_edge = find(diff(dummy) == 1) +1; % Location of first high value at rising edge
    falling_edge = find(diff(dummy) == -1);  % Location of last high value at falling edge
    
    epsilon = 0.1; 
    StartTime = (rising_edge/Fs) - epsilon;    % Convert to seconds
    StopTime = (falling_edge/Fs) + epsilon;

    View = "Spectrogram 1";
    View = repmat(View,height(StartTime),1);
    
    selection = (1:height(StartTime))';
    
    LowFreq = 900;
    LowFreq = repmat(LowFreq,height(StartTime),1);
    
    HighFreq = 3300;
    HighFreq = repmat(HighFreq,height(StartTime),1);
    
    Channel = ones(height(StartTime),1); % Only process channel 1
    
    Detected_Label = table(selection, View, Channel, StartTime, StopTime, LowFreq, HighFreq, ...
        VariableNames = ["Selection","View","Channel","Begin Time (s)","End Time (s)", "Low Freq (Hz)" ...
        "High Freq (Hz)"]); 

end

function Matched_Pred_Calls_Times = find_Suitable_Calls(Detect_Loc,Temp,t_avg,t_min,t_max,D_max) 
    % This function generate the search window for each detected calls
    v = 331.3+0.606*Temp;   % Speed of sound (meters/second)

    Matched_Pred_Calls_Times = cell(height(Detect_Loc),1);
    for i = 1:height(Detect_Loc)
        sz = [height(Detect_Loc{i,1}) 2];
        varTypes = ["double","double"];
        varNames = ["Start Time","End Time"];
        Matched_Pred_Calls_Times{i,1} = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
        for j = 1:height(Detect_Loc{i,1})
            W = [0,0];
            if t_min < (Detect_Loc{i,1}{j,5}-Detect_Loc{i,1}{j,4}) && (Detect_Loc{i,1}{j,5}-Detect_Loc{i,1}{j,4}) < t_avg
                W = [Detect_Loc{i,1}{j,4} - (t_avg - (Detect_Loc{i,1}{j,5}-Detect_Loc{i,1}{j,4}))/2 - D_max/v,...
                     Detect_Loc{i,1}{j,5} + (t_avg - (Detect_Loc{i,1}{j,5}-Detect_Loc{i,1}{j,4}))/2 + D_max/v];
    
            elseif t_avg <= (Detect_Loc{i,1}{j,5}-Detect_Loc{i,1}{j,4}) && (Detect_Loc{i,1}{j,5}-Detect_Loc{i,1}{j,4}) < t_max
                W = [Detect_Loc{i,1}{j,4} - D_max/v,...
                     Detect_Loc{i,1}{j,5} + D_max/v];
    
            end
            if W(1) ~= 0
                Matched_Pred_Calls_Times{i,1}{j,1} = W(1);
                Matched_Pred_Calls_Times{i,1}{j,2} = W(2);
            end
        end
    end
end

function Mask_Pred = generate_Label_Mask(Matched_Pred_Calls_Times,data_Test,Fs)
    % This function generate a search window mask for all detected calls
    % on each recorder
    Mask_Pred = cell(height(Matched_Pred_Calls_Times),1);

    for i = 1:height(Matched_Pred_Calls_Times)
        TestLabel = zeros(length(data_Test{1,1}),1); 
        
        for j = 1:height(Matched_Pred_Calls_Times{i,1})
            if Matched_Pred_Calls_Times{i,1}{j,1} ~= 0
                sl = Matched_Pred_Calls_Times{i,1}{j,1}; 
                sh = Matched_Pred_Calls_Times{i,1}{j,2};
                il = round(sl*Fs);
                ih = round(sh*Fs);
        
                TestLabel(il:ih,1) = 1; % For Train Data, 1Cha
            end
        end
        Mask_Pred{i,1} = TestLabel;
    end
end

function logical_Mask_Pred = generate_Logical_Label_Mask(Mask_Pred)
    % This function generate a logical search window mask for all 
    % detected calls on each recorder
    logical_Mask_Pred = Mask_Pred;
    for i = 1:height(Mask_Pred)
        dummy = Mask_Pred{i,1};
        logical_Mask_Pred{i,1} = logical(dummy); 
    end
end

function Matched_Pred_Calls = generate_Matched_Window(Intersect_Pred,Detect_Loc,data_Test)
    % This function find the set of clips containing matched calls
    sz = [height(Intersect_Pred) height(data_Test)+3];
    varTypes = ["cell","cell","cell","cell",...
        "cell","cell","cell","cell",...
        "cell","cell","cell","cell",...
        "double","double","double"];
    varNames = ["Mic 1 Call","Mic 2 Call","Mic 3 Call","Mic 4 Call",...
        "Mic 5 Call","Mic 6 Call","Mic 7 Call","Mic 8 Call",...
        "Mic 9 Call","Mic 10 Call","Mic 11 Call","Mic 12 Call",...
        "Start Time","End Time","Num Calls Assigned"];
    
    Matched_Pred_Calls = table('Size',sz,'VariableTypes',varTypes([1:height(data_Test),13:15]),'VariableNames',varNames([1:height(data_Test),13:15])); 
    
    for i = 1:height(Intersect_Pred)
        for j = 1:height(Detect_Loc)
            for k = 1:height(Detect_Loc{j,1})
                if (Detect_Loc{j,1}{k,4} >= Intersect_Pred{i,4}-1e-4) && (Detect_Loc{j,1}{k,5} <= Intersect_Pred{i,5}+1e-4)
                    Matched_Pred_Calls{i,j}{:,1} = [Matched_Pred_Calls{i,j}{:,1};k];
                    Matched_Pred_Calls{i,j} = Matched_Pred_Calls{i,j}(~cellfun('isempty',Matched_Pred_Calls{i,j}));
                    Matched_Pred_Calls{i,"Num Calls Assigned"} = Matched_Pred_Calls{i,"Num Calls Assigned"} + 1;
                    if Matched_Pred_Calls{i,"Start Time"} == 0
                        Matched_Pred_Calls{i,"Start Time"} = Detect_Loc{j,1}{k,4};
                    end
                    if Matched_Pred_Calls{i,"Start Time"} > Detect_Loc{j,1}{k,4}
                        Matched_Pred_Calls{i,"Start Time"} = Detect_Loc{j,1}{k,4};
                    end
                    if Matched_Pred_Calls{i,"End Time"} < Detect_Loc{j,1}{k,5}
                        Matched_Pred_Calls{i,"End Time"} = Detect_Loc{j,1}{k,5};
                    end
                end
            end     
        end
    end
end

function Matched_Pred_Calls = trim_Detection(Matched_Pred_Calls,data_Test)
    % Trim detection without 3+ detections across 4 locations
    TrimIndex = [];
    for i = 1:(height(Matched_Pred_Calls))
        not_zero_count = 0;
        for j = 1:height(data_Test)
            if ~isempty(Matched_Pred_Calls{i,j}{:})
                not_zero_count = not_zero_count + 1;
            end
        end
        if not_zero_count < 3
            TrimIndex = [TrimIndex;i];
        end    
    end
    Matched_Pred_Calls(TrimIndex,:) = [];
end

function Matched_Pred_Calls = seperate_Calls(Matched_Pred_Calls,Detect_Loc,t_avg,min_seperation_length,epsilon,data_Test,D_max,Temp)
    % Seperate > 2 matched calls grouped to 1 set of matched clips
    v = 331.3+0.606*Temp;   % Speed of sound (meters/second)
    varNames = ["Mic 1 Call","Mic 2 Call","Mic 3 Call","Mic 4 Call",...
        "Mic 5 Call","Mic 6 Call","Mic 7 Call","Mic 8 Call",...
        "Mic 9 Call","Mic 10 Call","Mic 11 Call","Mic 12 Call",...
        "Start Time","End Time","Num Calls Assigned"];    
    i = 1;
    while i <= height(Matched_Pred_Calls)
        if (Matched_Pred_Calls{i,"End Time"} - Matched_Pred_Calls{i,"Start Time"}) > 2*(t_avg+min_seperation_length)    % 2*(t_avg+D_max/v)
            multi_count = 0;
            for j = 1:height(Detect_Loc)
                if numel(Matched_Pred_Calls{i,j}{:}) > 1
                    multi_count = multi_count + 1;
                end
            end
            if multi_count >= 3
                % Find first call of the first grouped calls
                First_Group_ST = [];
                for j = 1:height(Detect_Loc)
                    if numel(Matched_Pred_Calls{i,j}{1}) >0
                        First_Group_ST = [First_Group_ST; j, Detect_Loc{j,1}{Matched_Pred_Calls{i,j}{1}(1),4}];
                    end
                end
                First_Detection_Index = find(First_Group_ST(:,2) == min(First_Group_ST(:,2)));
                Mic_Index = First_Group_ST(First_Detection_Index(1),1);
                Call_Index = Matched_Pred_Calls{i,Mic_Index}{1}(1);
                
                G_1 = cell(height(data_Test),1);
                G_2 = cell(height(data_Test),1);
                NumCalls_G_1 = 0;
                NumCalls_G_2 = 0;
                
                if Detect_Loc{Mic_Index,1}{Call_Index+1,4} - (Detect_Loc{Mic_Index,1}{Call_Index,5}+ min_seperation_length) > epsilon
    
                    W_1 = [Matched_Pred_Calls{i,"Start Time"}, Detect_Loc{Mic_Index,1}{Call_Index,5}+ min_seperation_length];
                    W_2 = [Detect_Loc{Mic_Index,1}{Call_Index,5}+ min_seperation_length+ epsilon, Matched_Pred_Calls{i,"End Time"}];
                    for j = 1:height(data_Test)
                        for k = 1:numel(Matched_Pred_Calls{i,j}{:})
                            Call_Index = Matched_Pred_Calls{i,j}{1}(k);
                            Mic_Index = j;
                            if Detect_Loc{Mic_Index,1}{Call_Index,4} < W_1(2)   % {Call_Index,4}
                                G_1{j,1} = [G_1{j,1},Call_Index];               % G_1{j,1} = Call_Index
                                NumCalls_G_1 = NumCalls_G_1 + 1;
                            end
                        end
    
                    end
                    for j = 1:height(Detect_Loc)
                        if numel(G_1{j,1}) > 0
                            dummy = Matched_Pred_Calls{i,j}{1}(:) ~= G_1{j,1};
                            index = all(dummy,2);
                            G_2{j,1} = Matched_Pred_Calls{i,j}{1}(index);
                            NumCalls_G_2 = NumCalls_G_2 + numel(G_2{j,1});
                        else
                            G_2{j,1} = Matched_Pred_Calls{i,j}{1}(:);
                            NumCalls_G_2 = NumCalls_G_2 + numel(G_2{j,1});
                        end
                    end
                    G_1_Row = cell2table([G_1',W_1(1),W_1(2),NumCalls_G_1],'VariableNames',varNames([1:height(data_Test),13:15]));
                    G_2_Row = cell2table([G_2',W_2(1),W_2(2),NumCalls_G_2],'VariableNames',varNames([1:height(data_Test),13:15]));
                    Matched_Pred_Calls = [Matched_Pred_Calls(1:(i-1),:);G_1_Row;G_2_Row;Matched_Pred_Calls((i+1):height(Matched_Pred_Calls),:) ];
                end
                
                
            end
        end
        i = i +1;
    end
end

function Matched_Pred_Vec = generate_Mask_Visualization(Matched_Pred_Calls,Intersect_Test_Vec,Fs)
    % Generate mask to visualize
    Matched_Pred_Vec = zeros(height(Intersect_Test_Vec),1);
    for i = 1:height(Matched_Pred_Calls)
    
        sl = Matched_Pred_Calls{i,"Start Time"}; %Train Data, 1 channel
        sh = Matched_Pred_Calls{i,"End Time"};
    
        il = round(sl*Fs);
        ih = round(sh*Fs);
    
        Matched_Pred_Vec(il:ih,1) = 1; % For Train Data, 1Cha
    end
end