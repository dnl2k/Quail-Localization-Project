
function [Detect_Loc, label_Pred] = Detection_Module(data_Test, Fs)
    % This function used pretrained DL model to detect the quail calls from
    % audio recordings
    %
    % Input:
    % 1. data_Test: the audio recordings
    % 2. Fs: sampling frequency of the inputted recordings
    % 
    % Output:
    % 1. Detect_Loc: table contain details of each detected calls,
    % formatted using Raven Pro table format to aid validation process.
    % 2. label_Pred: 1D vector contain the detection results for each
    % audio sample.

    label_Pred = cell(height(data_Test),1);

    % load Model
    [fPath, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
    cd([fPath,'\Model'])
    load('BiLSTM_Model_25-May-2023_23-07-30.mat','net');
    
    sequenceLength = 96e3;  
    sequenceOverlap = 0;   
    OVP = 0.5; 
    
    sequenceEdge = 1/4*sequenceLength;
    
    FrameSize = 800;
    OverlapPercent = 0.5;
    NumBand = 128;
    FrameOverlap = round(OverlapPercent*FrameSize);
    FFTLength = 2*FrameSize;
    LowF = 500;
    HighF = 4000;
    
    
    for k = 1:height(data_Test)
        
        % Preprocessing raw audio to neural network inputs   
        addpath([fPath,'\Functions'])       
        testFeatureCell = helperFeatureVector2Sequence(data_Test{k,1}',sequenceLength,sequenceOverlap);
                
        for i = 1:height(testFeatureCell)
            testFeatureCell{i,1} = melSpectrogram(testFeatureCell{i,1}',Fs, ...
                           'Window',hann(FrameSize,'periodic'), ...
                           'OverlapLength',FrameOverlap, ...
                           'FFTLength',FFTLength, ...
                           'NumBands',NumBand, ...
                           'FrequencyRange',[LowF,HighF]);
    
        
            testFeatureCell{i,1} = mag2db(testFeatureCell{i,1});
            mean_Feature = mean(testFeatureCell{i,1});
            std_Feature = std(testFeatureCell{i,1});
            testFeatureCell{i,1} = (testFeatureCell{i,1} -mean_Feature)./std_Feature;
            
            dummy = -1*ones(1,NumBand,1,width(testFeatureCell{i,1}));
            for j = 1:width(testFeatureCell{i,1})
                dummy(:,:,:,j) = reshape(testFeatureCell{i,1}(:,j),[1 NumBand 1]);
            end
            testFeatureCell{i,1} = dummy;
        end
           
        
        % Predict
        YPredicted = classify(net,testFeatureCell,'MiniBatchSize',16);
    
        % Reconstruct 
        % (Oversampling detection vector to the same size of the raw
        % audio wave)    
        ReLabel = -1*ones(1,height(data_Test{1,1}));
        il = 1;
        for i = 1:height(testFeatureCell)
            
            predSample = -1*ones(1,(FrameSize + FrameOverlap*(size(testFeatureCell{i,1},4)-1)));
            
            predSample(1,1:FrameSize) = double(YPredicted{i,1}(:,1));
        
            for j = 2:size(testFeatureCell{i,1},4)        
                li = ((j-1)*FrameOverlap+1);
                predSample(1,li:(li+FrameOverlap-1)) = double(YPredicted{i,1}(:,j));
            end
            ih = il + width(predSample)-1;
            
            predSample(1,1:sequenceEdge) = 0;
            predSample(1,(width(predSample)-sequenceEdge+1):width(predSample)) = 0;
    
            ReLabel(1,il:ih) = predSample;
            il = ih +1;
        end
        label_Pred{k,1} = ReLabel;
        
        %------ Make overlapped Test Data Sequence --------% 
        % (Since the recordings were first seperated into non-overlaping
        % 4-second clips then the network detect calls on those clips, 
        % it is possible that a covey call is located between 2 consecutive
        % clips, e.g. half of a call on the n-th clip and the remaining half
        % on the n+1-th clip. This overlapped detection step is to mitigate
        % the above issue.)

        % Preprocessing raw audio to neural network inputs
        Num_OVP = OVP*sequenceLength;
        testFeatureCell = helperFeatureVector2Sequence(data_Test{k,1}((Num_OVP+1):height(data_Test{k,1}),:)',sequenceLength,sequenceOverlap);      
        for i = 1:height(testFeatureCell)
            testFeatureCell{i,1} = melSpectrogram(testFeatureCell{i,1}',Fs, ...
                           'Window',hann(FrameSize,'periodic'), ...
                           'OverlapLength',FrameOverlap, ...
                           'FFTLength',FFTLength, ...
                           'NumBands',NumBand, ...
                           'FrequencyRange',[LowF,HighF]);
    
        
            testFeatureCell{i,1} = mag2db(testFeatureCell{i,1});
            mean_Feature = mean(testFeatureCell{i,1});
            std_Feature = std(testFeatureCell{i,1});
            testFeatureCell{i,1} = (testFeatureCell{i,1} -mean_Feature)./std_Feature;
            
            dummy = -1*ones(1,NumBand,1,width(testFeatureCell{i,1}));
            for j = 1:width(testFeatureCell{i,1})
                dummy(:,:,:,j) = reshape(testFeatureCell{i,1}(:,j),[1 NumBand 1]);
            end
            testFeatureCell{i,1} = dummy;
        end
           
        % Predict overlapped
        YPredicted = classify(net,testFeatureCell,'MiniBatchSize',16);
        
        % Reconstruct    
        ReLabel = -1*ones(1,height(data_Test{1,1}));
        il = Num_OVP+1;
        for i = 1:height(testFeatureCell)
            
            predSample = -1*ones(1,(FrameSize + FrameOverlap*(size(testFeatureCell{i,1},4)-1)));
            
            predSample(1,1:FrameSize) = double(YPredicted{i,1}(:,1));
        
            for j = 2:size(testFeatureCell{i,1},4)        
                li = ((j-1)*FrameOverlap+1);
                predSample(1,li:(li+FrameOverlap-1)) = double(YPredicted{i,1}(:,j));
            end
            ih = il + width(predSample)-1;
            
            predSample(1,1:sequenceEdge) = 0;
            predSample(1,(width(predSample)-sequenceEdge+1):width(predSample)) = 0;
    
            ReLabel(1,il:ih) = predSample;
            il = ih +1;
        end
    
        % Combine the results from overlapped
        label_Pred{k,1} = label_Pred{k,1} + ReLabel(:,1:width(label_Pred{k,1}));
        index = label_Pred{k,1}(1,:) >= 2;
        label_Pred{k,1}(:,~index) = 1;
    end
    
    % Write Predicted Table in Raven Pro format 
    Detect_Loc = cell(height(label_Pred),1);
    
    for i = 1:height(label_Pred)
        
        dummy = label_Pred{i,1};
        index = dummy(1,:) == 2; % Quail A (covey call)         
    
        dummy(:,index) = 2;
        dummy(:,~index) = 1;
    
        rising_edge = find(diff(dummy) == 1) +1; % Location of first high value at rising edge
        falling_edge = find(diff(dummy) == -1);  % Location of last high value at falling edge
        
        StartTime = (rising_edge/Fs)';    % Convert to seconds
        StopTime = (falling_edge/Fs)';
    
        View = "Spectrogram 1";
        View = repmat(View,height(StartTime),1);
    
        selection = (1:height(StartTime))';
    
        LowFreq = 900;
        LowFreq = repmat(LowFreq,height(StartTime),1);
    
        HighFreq = 3300;
        HighFreq = repmat(HighFreq,height(StartTime),1);
    
        Channel = ones(height(StartTime),1); % Only process channel 1
        
        Detect_Loc{i,1} = table(selection, View, Channel, StartTime, StopTime, LowFreq, HighFreq, ...
            VariableNames = ["Selection","View","Channel","Begin Time (s)","End Time (s)", "Low Freq (Hz)" ...
            "High Freq (Hz)"]); 
    end
end
