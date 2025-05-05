%% Section 1. Input for script to analyze movies of ARF condensates

% Using the TrackMate trajectories, fit an ellipse to each condensate

% Output of the second section is a .csv file with the time series of the
% fits to each condensate

% The third section is for plotting

videoname = 'Individual_1.avi';
time_resolution = 0.076; % [s]
pixel_resolution = 0.304; % [µm]
window = 5; % Fixed window in x and y around center pixel. Depends on spot size.

linkdata = readmatrix('TrackMate_data_for_Individual_1.xlsx','Sheet','Links');
% Second column is track IDs, sixth column is edge times, seventh and eighth
% columns are x and y data

%% Section 2. Construct a time series of the fits to each condensate

video = VideoReader(videoname);

possible_edge_times = (1:video.NumFrames).*time_resolution - time_resolution/2;
% This analysis ignores TrackMate data that doesn't correspond to a frame at a known time point.

% Extract track IDs, edge times, and x and y data
[unique_entries, ~, indices] = unique(linkdata(:,2), 'stable');
% Ensure proper alignment
times = cell(size(unique_entries));
x_positions = cell(size(unique_entries));
y_positions = cell(size(unique_entries));
for i = 1:length(times)
    times{i} = NaN;
    x_positions{i} = NaN;
    y_positions{i} = NaN;
end
% Iterate through the unique entries and assign the corresponding times
for i = 1:length(unique_entries)
    idx = find(indices == i);
    times{i} = linkdata(idx, 6);
    x_positions{i} = linkdata(idx, 7);
    y_positions{i} = linkdata(idx, 8);    
end

% Check if a condensate is identified in a frame
ncondensates = length(unique(linkdata(:,2)));
check_if_in_frame = zeros(length(possible_edge_times),ncondensates);
frame = cell(length(ncondensates),1);
tol = 1e-3; % Tolerance (used later, too)
for ii = 1:ncondensates
    check_if_in_frame(:,ii) = ismembertol(possible_edge_times,times{ii},tol,'DataScale',1);
    frame{ii} = find(check_if_in_frame(:,ii));
end

ii = 0;
idx2 = zeros(video.NumFrames,ncondensates);
w1 = cell(video.NumFrames,ncondensates);
x_frame = zeros(video.NumFrames,ncondensates);
y_frame = zeros(video.NumFrames,ncondensates);
x_for_ellipse = cell(video.NumFrames,ncondensates);
y_for_ellipse = cell(video.NumFrames,ncondensates);
Cent = cell(video.NumFrames,ncondensates);
MajAx = zeros(video.NumFrames,ncondensates);
MinAx = zeros(video.NumFrames,ncondensates);
Orient = zeros(video.NumFrames,ncondensates);
channel_to_read = zeros(video.Height,video.Height,video.NumFrames);

% Analyze all frames
while hasFrame(video)

    ii = ii+1;
    frame_to_read = readFrame(video); % Read all frames
    channel_to_read(:,:,ii) = frame_to_read(:,:,3); % Channel 3 for ARFs
    image_to_analyze = channel_to_read(:,:,ii);

    for jj = 1:ncondensates
        
        check_condensate = check_if_in_frame(:,jj);

        if check_condensate(ii) == 1

            edge_x = x_positions{jj};
            edge_y = y_positions{jj};
            edge_times_to_check = times{jj};
            idx1 = find(abs(edge_times_to_check-possible_edge_times(ii))<tol);

            x_center_orig = round(edge_x(idx1)/pixel_resolution);
            y_center_orig = round(edge_y(idx1)/pixel_resolution);

            x_frame(ii,jj) = x_center_orig;
            y_frame(ii,jj) = y_center_orig;
    
            % Do not analyze if the condensate is too close to the edge
            if x_center_orig > window && y_center_orig > window
                if x_center_orig < (video.Height - window) && y_center_orig < (video.Height - window)
                    
                    % Keep track of indices
                    ind_to_keep = idx1;
                    idx2(ii,jj) = ind_to_keep;

                    % Localize
                    w = image_to_analyze(y_center_orig-window:y_center_orig+window,...
                        x_center_orig-window:x_center_orig+window); % Note order based on TrackMate data

                    w1{ii,jj} = w; % Will need for plotting

                    % Tighten window at the periphery of the image since
                    % those data are presumably from another condensate
                    w = padarray(w(2:end-1,2:end-1),[1 1]);
                    % (Equivalently, we could just use a smaller window at
                    % the start, but this is better for plotting.)

                    % Binarize image via global thresholding of 33% of max
                    BW = imbinarize(w,max(max(w))/3);

                    % In a few cases, a nearby condensate will also appear
                    % in the window. Assume the condensate-to-analyze is
                    % nearest to the center pixel.
                    CC = bwconncomp(BW);
                    C = CC.PixelIdxList;
                    [nx,ny] = find(BW);
                    distance = sqrt((window+1-nx).^2+(window+1-ny).^2);
                    A1 = cat(1,C{:});
                    pixel_value = A1(distance == min(distance));
                    rem = cellfun(@(x) find(x == pixel_value(1)),C,'un',0);
                    cell_idx = find(~cellfun('isempty',rem));
                    A2 = C{cell_idx};
                    center_component = zeros(2*window+1,2*window+1);
                    center_component(A2) = 1;

                    % Calculate centroid, orientation and major/minor axis length of the ellipse
                    s = regionprops(center_component,{'Centroid','Orientation','MajorAxisLength','MinorAxisLength'});
            
                    % Calculate the ellipse "line"
                    theta = linspace(0,2*pi);
                    col = (s.MajorAxisLength/2)*cos(theta);
                    row = (s.MinorAxisLength/2)*sin(theta);
                    M = makehgtform('translate',[s.Centroid, 0],'zrotate',deg2rad(-1*s.Orientation));
                    D = M*[col;row;zeros(1,numel(row));ones(1,numel(row))];
                    x_for_ellipse{ii,jj} = D(1,:);
                    y_for_ellipse{ii,jj} = D(2,:);
            
                    % Statistics
                    Cent{ii,jj} = s.Centroid;
                    MajAx(ii,jj) = s.MajorAxisLength;
                    MinAx(ii,jj) = s.MinorAxisLength;
                    Orient(ii,jj) = s.Orientation;

                    if s.Orientation == 0
                        Orient(ii,jj) = 0.01; % Dummy number to be replaced later
                    end

                end

            end

        end

    end
        
end

% Output data from fit as a function of time

edge_times_all = cell(ncondensates,1);
majAx_all = cell(ncondensates,1);
minAx_all = cell(ncondensates,1);
orient_all = cell(ncondensates,1);

for jj = 1:ncondensates

    edge_times_to_check = times{jj};
    edge_times_to_use = edge_times_to_check(nonzeros(idx2(:,jj)));
    majAx_to_keep = nonzeros(MajAx(:,jj));
    minAx_to_keep = nonzeros(MinAx(:,jj));
    orient_to_keep = nonzeros(Orient(:,jj));

    % Add NaN's wherever a frame was skipped
    dum_var = ismember(times{jj},edge_times_to_use).*times{jj};
    dum_var(dum_var==0) = NaN;
    edge_times_all{jj} = dum_var;
    dum_var1 = dum_var;
    dum_var2 = dum_var;
    dum_var1(dum_var1>0) = majAx_to_keep;
    majAx_all{jj} = dum_var1;
    dum_var(dum_var>0) = minAx_to_keep;
    minAx_all{jj} = dum_var;
    dum_var2(dum_var2>0) = orient_to_keep;
    dum_var2(dum_var2==0.01) = 0; % Replace earlier dummy number in orientations with zeros
    orient_all{jj} = dum_var2;

end

% Write to .csv file

edge_times = cell2mat(edge_times_all);
majorAxisLength = cell2mat(majAx_all);
minorAxisLength = cell2mat(minAx_all);
orientationOfEllipse = cell2mat(orient_all);

% Condensate ID, times, x and y positions, major and minor axes, aspect
% ratio, area, orientation of ellipse
A_out = [linkdata(:,2),edge_times,linkdata(:,7),linkdata(:,8),majorAxisLength*pixel_resolution,...
    minorAxisLength*pixel_resolution,majorAxisLength./minorAxisLength,...
    pi*(majorAxisLength/2).*(minorAxisLength/2).*pixel_resolution^2,orientationOfEllipse];
T = array2table(A_out);
T.Properties.VariableNames(1:size((A_out),2)) = {'Condensate_ID','Times (s)','x_position (um)',...
    'y_position (um)','Length_of_major_axis (um)','Length_of_minor_axis (um)','Aspect_ratio',...
    'Area (um^2)','Angle (degrees)'};
writetable(T,'image_analysis_output.csv')

%% Section 3. Plotting

% Plot an entire frame

frame_to_plot = 730; % Change this depending on what's in the frame

figure(1);
imagesc(channel_to_read(:,:,frame_to_plot))
cmax = round(max(max(channel_to_read(:,:,frame_to_plot))));
colorbar
colormap('hot')
caxis([0 cmax])
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'fontsize',24)

hold off

% Identify condensates

N_ID = find(idx2(frame_to_plot,:)); % Condensate IDs

figure(2);
imagesc(channel_to_read(:,:,frame_to_plot))
colorbar
colormap('hot')
caxis([0 cmax])
hold on
x_to_plot = x_frame(frame_to_plot,N_ID);
y_to_plot = y_frame(frame_to_plot,N_ID);
for kk = 1:length(N_ID)
    rectangle('Position',[x_to_plot(kk)-window y_to_plot(kk)-window 2*window+1 2*window+1],'EdgeColor','w','LineWidth',2);
end
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'fontsize',24)

hold off

% Plot zoomed-in images of all condensates in the frame

for kk = 1:length(N_ID)
    figure(2+kk);
    imagesc(w1{frame_to_plot,N_ID(kk)})
    colorbar
    colormap('hot')
    caxis([0 cmax])
    hold on
    plot(x_for_ellipse{frame_to_plot,N_ID(kk)},y_for_ellipse{frame_to_plot,N_ID(kk)},'w','LineWidth',4)
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    set(gca,'fontsize',24)
    hold off
end

hold off
