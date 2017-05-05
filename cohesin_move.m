function cohesin_move(infile,outfile)
%This program will parse an outfile identify an unlinked cohesin ring and
%move place it around the nearest DNA bead.  This program was written for
%cohesin rings that become displaced after condensin tethering of distal
%DNA beads causes cohesin rings to fall off the strand.

%% Load in the .out file
%get fileID
fid = fopen(infile);
tline = fgetl(fid);
%loop through the file line by line
%with this setup the end-of-file (EOF) is -1, so test for character to ID
%EOF
%instantiate color_tog switch to 0
color_tog = 0;
%instantiate time_count
time_count = 0;
%instantiate time_tog
time_tog = 0;
while ischar(tline)
    %% Parse .out file for color list
    %if statment to turn off toggle if it finds a blank line
    if strcmp(tline,'') == 1 && color_tog == 1
        color_tog = 0;
        time_tog = 1;
    end
    %if statement to parse colors
    if color_tog == 1
        colors(color_count,1) = str2double(tline);
        color_count = color_count + 1;
    end
    %if statement to turn on color_tog if MassColor line is found
    if strcmp(tline,'MassColors') == 1
        color_tog = 1;
        color_count = 1;
    end
    %% Parse the time coords
    %if statement that will parse coordinates
    %time_tog must be on (indicating color section is finished)
    %line must not be blank, blank lines occurs between time points in file
    %line must not be a "Time" line, only lines with 2 entries are "Time"
    %lines
    if time_tog == 1 && strcmp(tline,'') == 0 && length(strsplit(tline)) ~= 2
        coords = strsplit(tline);
        coords = [str2double(coords{1}),...
            str2double(coords{2}),...
            str2double(coords{3})];
        coord_mat(coord_count,:,time_count) = coords;
        coord_count = coord_count + 1;
    end
    %if statment beginning time_count variable and coord_count variable
    %checks if tline is a "Time" line and that color section is finished
    if length(strsplit(tline)) == 2 && time_tog == 1
        time_count = time_count + 1;
        coord_count = 1;
        time_line{time_count} = tline;
    end
    
    
    %get next line of file
    tline = fgetl(fid);
end
%close fid
fclose(fid);
%% Calculate the mean of position of each XYZ position
%cohesin color number is 5, so create binary array with 5 = 1
coh_bin = colors == 5;
%apply binary to first and last timepoint
coh_mat_final = [coord_mat(coh_bin,1,end),...
    coord_mat(coh_bin,2,end),...
    coord_mat(coh_bin,3,end)];
coh_mat_init = [coord_mat(coh_bin,1,1),...
    coord_mat(coh_bin,2,1),...
    coord_mat(coh_bin,3,1)];
%% Find min distance from DNA beads for each ring
%divide coh_mat_final to each ring (each ring has 16 beads)
ring_num = size(coh_mat_final,1)/16;
%loop through each ring coords in X,Y, and Z and calc mean coordinate
for i = 1:ring_num
    ring_mat_final(i,:) = [mean(coh_mat_final(1+(16*(i-1)):16+(16*(i-1)),1))...
        mean(coh_mat_final(1+(16*(i-1)):16+(16*(i-1)),2)),...
        mean(coh_mat_final(1+(16*(i-1)):16+(16*(i-1)),3))];
    ring_mat_init(i,:) = [mean(coh_mat_init(1+(16*(i-1)):16+(16*(i-1)),1))...
        mean(coh_mat_init(1+(16*(i-1)):16+(16*(i-1)),2)),...
        mean(coh_mat_init(1+(16*(i-1)):16+(16*(i-1)),3))];
end
%create binary for dna coordinates
dna_bin = ~coh_bin;
%parse the final timepoint for all dna bead positions
dna_mat_final = [coord_mat(dna_bin,1,end),...
    coord_mat(dna_bin,2,end),...
    coord_mat(dna_bin,3,end)];
dna_mat_init = [coord_mat(dna_bin,1,1),...
    coord_mat(dna_bin,2,1),...
    coord_mat(dna_bin,3,1)];
%calculate the Euclidean distance between each ring and each dna
%bead
for j = 1:size(dna_mat_final,1)
    for k = 1:ring_num
        dist_mat_final(j,k) = sqrt(sum((dna_mat_final(j,:) - ring_mat_final(k,:)).^2));
        dist_mat_init(j,k) = sqrt(sum((dna_mat_init(j,:) - ring_mat_init(k,:)).^2));
    end
end
%find the minimum distance for each ring
min_dist_ring = min(dist_mat_final);
%the largest of the min distances is the displaced ring
%take the max of mins to find ring idx
[~,ring_idx] = max(min_dist_ring);

%% Identify which dna bead was closest to displaced cohesin ring at start
disp_ring_dist = dist_mat_init(:,ring_idx);
[~,dna_idx] = min(disp_ring_dist);
%parse that dna bead's current coordinates
prox_bead_coords = dna_mat_final(dna_idx,:);

%% move cohesin ring back around prox. dna bead coordinate set
%parse the displaced ring's coordinates
ring_total_idx = 1+((ring_idx-1)*16):16*(ring_idx);
ring_coords = coh_mat_final(ring_total_idx,:);
%Find the indices of the correct coordinates
coh_idxs = find(coh_bin);
%This poriton of code is redundant if cohesin rings are first in color list
disp_coh_idxs = coh_idxs(ring_total_idx);
for m = 1:size(ring_coords,1)
    coord_mat(disp_coh_idxs(m),:,end) = ring_coords(m,:) - ring_mat_final(ring_idx,:)...
        + prox_bead_coords;
end

%% Write out new .out file with moved cohesin
%open file to write
fid_out = fopen(outfile,'w');
%re-open the old .out file
fid_in = fopen(infile);
%grab the first line of old .out file
tline = fgetl(fid_in);
%loop through
while ischar(tline)
    %when first timeline is reached
    if strcmp(tline,time_line{1}) == 1
        %print out the final timeline
        fprintf(fid_out,'%s\n',time_line{end});
        %print out the altered coordinates
        for k = 1:size(coord_mat,1)
            fprintf(fid_out,'%g %g %g\n',...
                coord_mat(k,1,end),coord_mat(k,2,end),coord_mat(k,3,end));
        end
        %Set end-of-line after printing all coords
        tline = -1;
    else
        fprintf(fid_out,'%s\n',tline);
        tline = fgetl(fid_in);
    end
end
fclose(fid_in);
fclose(fid_out);
end
