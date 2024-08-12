% Return the cumulator matrix and range of the matrix 
% Improve scalability
% Fix bugswhen increasing the number of mic/tau inputs

function [cum_result_matrix,range,localization_log, Pred_Loc] = Localization_Module(tau, tau_corr, Mic, Temp, resolution, range_offset)
    % This function performs hyperbolic localization using algorithm
    % from Quail Localization paper.
    % Input:
    % 1. tau: TDOA values
    % 2. tau_corr: Correlation coefficients corresponding to each TDOA
    % values
    % 3. Mic: Recorders' location
    % 4. Temp: Air temperature, in degree Celcius.
    % 5. resolution: the resolution of the accumulator matrix, in meters
    % 6. range_offset: extensions from the corners of array, in meters 
    % (assuming square array, will need to adjust the code for different 
    % array's geometry)
    % 
    % Output:
    % 1. cum_result_matrix: the result accumulator matrix
    % 2. range: the range of the accumulator matrix, used to support
    % plotting the cum_result_matrix if needed
    % 3. localization_log: log to report about whether the localization was
    % success or not for each matched calls.
    % 4. Pred_Loc: list of the coordinates of the localized calls.
    
    sz = [height(tau) 3];
    varTypes = ["double","string","string"];
    varNames = ["Time lag#","Status","Error Note"];
    localization_log = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
    localization_log{:,1} = (1:height(tau))';

    c=331.3+0.606*Temp;     % Speed of sound formula

    tau_corr = abs(tau_corr);   % Making sure the corr is positive

    range=[min(Mic)-range_offset;max(Mic)+range_offset];    % Set up the range for the meshgrid

    [X, Y] = meshgrid(range(1,1):resolution:range(2,1), range(1,2):resolution:range(2,2));  % Initialize the meshgrid
    cum_result_matrix = zeros(size(X));     % Initialize the final result matrix for all matched calls
    Pred_Loc = [];                          % To store Predicted Location

    for i = 1:height(tau)
        result_matrix = zeros(size(X));     % Initial the result matrix for each matched call
        Num_Valid = sum(double(~isnan(tau(i,:))));  % Find the valid tau values
        tau(i,:) = tau(i,:) - min(tau(i,:));    % Making sure the tau values are >= 0
        [tau_sorted,j]=sort(tau(i,:));          % Sort the tau, Mic loc., corr based on the tau values
        Mics=Mic(j,:);
        tau_corr_sorted = tau_corr(i,j);
        for m=1:Num_Valid-1                 % Filling the matrix with the hyperbolas
            for n = m+1:Num_Valid
                
                populated_matrix = (sqrt((X-Mics(m,1)).^2+(Y-Mics(m,2)).^2)- ...
                    sqrt((X-Mics(n,1)).^2+(Y-Mics(n,2)).^2))-c*(tau_sorted(m)-tau_sorted(n));                               

                M = contourc(populated_matrix,[0 0]);
                index = M(1,:) == 0;
                M(:,index) = [];
                index = sub2ind(size(populated_matrix),round(M(2,:)) ,round(M(1,:)));
                populated_matrix(:) = 0;
                populated_matrix(index) = 1;
                
                populated_matrix = imdilate(populated_matrix,strel("disk",1));
        
                result_matrix = result_matrix + populated_matrix;  
            end
        end
        
        
        %--------Binarize the accumulated image-----------% 
        max_result_matrix = max(result_matrix,[],'all');
        if max_result_matrix >= ceil((Num_Valid+1)/2)   % When at least the majority of TOA values is correct
            index = result_matrix >= ceil((Num_Valid+1)/2);
            result_matrix(index) = 1;
            result_matrix(~index) = 0;
        else                                            
            localization_log{i,"Status"} = "Failed";
            localization_log{i,"Error Note"} = "Invalid lag values";
            continue
        end
        %-------------------------------------------------%
        
        blobs = bwconncomp(result_matrix);  % Finding the blob
        if blobs.NumObjects < 1     % Checking if any blob exist
            localization_log{i,"Status"} = "Failed";
            localization_log{i,"Error Note"} = "No blob found";
            continue
        end
        stats = regionprops(blobs,'Centroid');
        centroid = vertcat(stats.Centroid);
        index = sub2ind(size(X),round(centroid(:,2)),round(centroid(:,1)));
        centroid = [X(index) Y(index)];
        
    
        if height(centroid) > 1
            TDOA_table = -1*ones(height(centroid), Num_Valid);
            for j = 1:height(centroid)
                dist = sqrt(sum((Mics(1:Num_Valid,:) - centroid(j,:)).^2,2));
                travel_time = dist./c;
                TDOA = travel_time - min(travel_time);
                TDOA_table(j,:) = TDOA';
            end             
            tau_diff = TDOA_table - tau_sorted(1,1:Num_Valid);
            
            dummy = tau_sorted(1,1:Num_Valid);
            index = dummy == 0;
            dummy(index) = 1e-100;
            % Scaled error with scale = original tau values
            tau_error = sum(abs(tau_diff(:,1:Num_Valid)./dummy.*tau_corr_sorted(1,1:Num_Valid)),2);

            for j = 1:height(tau_error)     % Region checking
                if TDOA_table(j,1) ~= 0
                    tau_error(j,1) = nan;
                end
            end

            index = find(tau_error == min(tau_error));
            Pred_Loc = [Pred_Loc;centroid(index,:)];
            index_blob = blobs.PixelIdxList{1,index};
            result_matrix(:) = 0;
            result_matrix(index_blob) = 1;
        else
            Pred_Loc = [Pred_Loc;centroid];
        end  
        
        localization_log{i,"Status"} = "Success";
        localization_log{i,"Error Note"} = "NA";
        cum_result_matrix = cum_result_matrix + result_matrix;

    end

end