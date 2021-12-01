%% Read ASL LIST/DATA
% Ananth Madhuranthakam
% Original - June 18, 2012
% Modified - Jan. 17, 2013
% Modified by Yue - October 28, 2013: replace zip_factor*nx to nx_recon2, zip_factor*ny_recon to ny_recon2
% Go to Line 117 to change the values based on reconstruction matrix in .list file.
% INPUT: 
% OUTPUT:
%% Read DATA file and setup parameters for recon
clear all;
close all;

tic
current_dir='../MATLAB/Data/2D Data/20201202_Kidney_ASL';
zip_factor = 2;
debug_mode = 1;
M0_series_no = 12;
use_M0 = 0;
write_data = 0;
write_quant_data = 1;
slice_thickness = 15.0;
some_thing_funny = 0;
sagittal = 0;       % 1 for sagittal and 0 for coronal

if (~use_M0)
    write_quant_data = 0;
end

[fname, dir_name] = uigetfile('*.data', 'Choose PCASL DATA file',current_dir);
filename = [dir_name fname];
[rawdata, info, ky_range] = readListData_YRV3(filename);
% [rawdata, info, ky_range] = readListData_YW_2(filename,0.7,0.5); %LM 2020-02-05
series_no = str2double(fname(6:7));

if  (ky_range(1) + ky_range(2)) > 0
    full_kspace = false;
else
    full_kspace = true;
end

%% Distribute all parameters
order = {'kx','ky','kz','loca','dyn','card','echo','mix','aver','extr1','chan'};
size_data = size(rawdata);
length_datasize = length(size_data);

nx_oversample = 2;
nx      = size_data(1)/nx_oversample;
if (length_datasize > 1)
    ny      = size_data(2);
else
    ny = 1;
end
if (length_datasize > 2)
    nz      = size_data(3);
else
    nz = 1;
end
if (length_datasize > 3)
    nlocs   = size_data(4);
else
    nlocs = 1;
end
if (length_datasize > 4)
    ndynamics  = size_data(5);
else
    ndynamics = 1;
end
if (length_datasize > 5)
    ncards  = size_data(6);
else
    ncards = 1;
end
if (length_datasize > 6)
    nechoes = size_data(7);
else
    nechoes = 1;
end
if (length_datasize > 7)
    nmix    = size_data(8);
else
    nmix = 1;
end
if (length_datasize > 8)
    navgs   = size_data(9);
else
    navgs = 1;
end
if (length_datasize > 9)
    lbl_ctrl= size_data(10);
else
    lbl_ctrl = 1;
end
if (length_datasize > 10)
    ncoils  = size_data(11);
else
    ncoils = 1;
end

%% Read M0 data for 2D scans
if (use_M0)
    [M0_fname, M0_dir_name] = uigetfile('*.data', 'Choose M0 DATA file',current_dir); % Yue adds
    M0_filename = [M0_dir_name M0_fname]; % Yue adds
    [M0_rawdata, M0_info, M0_ky_range] = readListData_YRV3(M0_filename); % Yue mods M0_filename
    ny_recon_M0 = 4*ceil((max(abs(M0_ky_range)))/2);
    size_M0_data = size(squeeze(M0_rawdata));
    nx_M0= size_M0_data(1)/nx_oversample;
end

% Determine ny_recon from ky_range
ny_recon = 4*ceil((max(abs(ky_range)))/2);

% ny_recon=160;
% ny_recon2 = 224; % Yue: round(FOV_PE*RFOV/recon_voxel_size) 360*113.3333%/2.32 (2.32 from ACQ voxel MPS in ExamCard Info)
% nx_recon2 = 240; % Yue: round(FOV_RO/recon_voxel_size); 360/2.31 from ACQ voxel MPS
% zip_factor_ny = ny_recon2/ny_recon2; % Yue


%% Find ny_recon2, nx_recon2 - Josh
list_name = strcat(filename(1:length(filename)-4),'list');
fileID = fopen(list_name,'r');

C = textscan(fileID, '%s', 'Delimiter', '\n');

y = strfind(C{1}, 'Y-resolution');
x = strfind(C{1}, 'X-resolution');
yrow = find(~cellfun('isempty', y));
xrow = find(~cellfun('isempty', x));
ystr=C{1}(yrow);
xstr=C{1}(xrow);

ny_recon2 = str2num(ystr{1}(length(xstr{1})-3:end));
nx_recon2 = str2num(xstr{1}(length(xstr{1})-3:end));
zip_factor_ny = ny_recon2/ny_recon; % Yue


fclose(fileID);
%%
% Squeeze data and the end result is:
% 2D: kx-ky-ndynamics-lbl_ctrl-ncoils (Need to modify for >1 slice)
% 3D: kx-ky-kz-lbl_ctrl-ncoils
data = squeeze(rawdata);
% sel_imgs=[1     6    11    13    14];
% sel_imgs=[6 13 14];
% data = data(:,:,sel_imgs,:,:);
% data = data(:,:,:,:,:);

if (use_M0)
    % M0 2D: kx-ky-lbl_ctrl-ncoils -> kx-ky-ncoils
    M0_data = squeeze(M0_rawdata);
    
    ny_M0      = size(M0_data,2);
%     size_data(2);
    zip_factor_ny_M0 = ny_recon2/ny_recon_M0; % Yue

    % M0_data = squeeze(mean(M0_data, 3));
    %     M0_data = squeeze(M0_data(:,:,1,:)); %% This was added to correct size for specific data set. Remove for subject 1 - Josh
end



%%
% Average control data and label data for each coil
ctrl_data = zeros(nx*nx_oversample,ny,nz,ncoils);
lbl_data = zeros(nx*nx_oversample,ny,nz,ncoils);

if (nz > 1) % 3D
    ctrl_data = squeeze(data(:,:,:,1,:));
    lbl_data = squeeze(data(:,:,:,2,:));
    nz_zip_factor = zip_factor;
else
    temp_data = squeeze(mean(data(:,:,:,1,:),3));
    ctrl_data = permute(temp_data, [1 2 4 3]);
    
    temp_data = squeeze(mean(data(:,:,:,2,:),3));
    lbl_data = permute(temp_data, [1 2 4 3]);  % average 4 dyns
       
%     temp_data = squeeze(data(:,:,1,:));
%     ctrl_data = permute(temp_data, [1 2 4 3]);
%     
%     temp_data = squeeze(data(:,:,2,:));
%     lbl_data = permute(temp_data, [1 2 4 3]);  %  LZ
    
    nz_zip_factor = 1;
end

sub_data = ctrl_data - lbl_data;

% Multiply with Fermi filter
% % Sinc filter
% len1 = 30;
% filt_cos = cos(pi*(linspace(10,0,len1+1)/20)).^2;
% low_filterX_1d = ones(1,nx*nx_oversample);
% low_filterX_1d(1:1+len1) = filt_cos;
% low_filterX_1d(end-len1:end) = filt_cos(end:-1:1);
%
% low_filterY_1d = ones(1,ny);
% low_filterY_1d(end-len1:end) = filt_cos(end:-1:1);
%
% low_filter_2d = low_filterX_1d'*low_filterY_1d;

if ~full_kspace % if partial kspace 
    % Fermi filter
    %YR START
    WindA = 16;
    WindB = 16;
    Ra = nx*nx_oversample/2;
    Rb = ny_recon/2;
    %YR END
    fermi_filter = fermi2d(nx*nx_oversample,ny_recon,Ra,Rb,WindA,WindB);
    low_filter_2d = fermi_filter(:,ny_recon-ny+1:end);
    if (use_M0)
        fermi_filter_M0 = fermi2d(nx_M0*nx_oversample,ny_recon_M0,nx_M0*nx_oversample/2,ny_recon_M0/2,16,16);
        low_filter_2d_M0 = fermi_filter_M0(:,ny_recon_M0-ny_M0+1:end);
    end 
    for (coil_ct = 1:ncoils)
        for (slice_ct = 1:nz)
            ctrl_data(:,:,slice_ct,coil_ct) = ctrl_data(:,:,slice_ct,coil_ct).*low_filter_2d;
            sub_data(:,:,slice_ct,coil_ct) = sub_data(:,:,slice_ct,coil_ct).*low_filter_2d;
%             ctrl_data(:,:,slice_ct,coil_ct) = ctrl_data(:,:,slice_ct,coil_ct);
%             sub_data(:,:,slice_ct,coil_ct) = sub_data(:,:,slice_ct,coil_ct);
            if (debug_mode)
                lbl_data(:,:,slice_ct,coil_ct) = lbl_data(:,:,slice_ct,coil_ct).*low_filter_2d;
%                lbl_data(:,:,slice_ct,coil_ct) = lbl_data(:,:,slice_ct,coil_ct);
            end
        end
        if (use_M0)
            % M0_data = .....*low_filter_2D_M0;
            M0_data(:,:,coil_ct) = M0_data(:,:,coil_ct).*low_filter_2d_M0;
%             M0_data(:,:,coil_ct) = M0_data(:,:,coil_ct);
        end
    end
    
    % Create Chopping matrix
    chopx = ones(nx_recon2, ny_recon2);
    chopy = ones(nx_recon2, ny_recon2);
    chopx(:,2:2:end) = 1*chopx(:,2:2:end);
    chopy(2:2:end,:) = 1*chopy(2:2:end,:);
    
    % img = zeros(nx_recon2, ny_recon2, nz*nz_zip_factor, ncoils);
    % low_res_img = zeros(nx_recon2, ny_recon2, nz*nz_zip_factor, ncoils);
    % if (use_M0)
    %     M0_img = zeros(nx_recon2, ny_recon2, ncoils);
    % end
    
    % For debugging purpose
    if (debug_mode)
        ctrl_img = zeros(nx_recon2, ny_recon2, nz*nz_zip_factor, ncoils);
        lbl_img = zeros(nx_recon2, ny_recon2, nz*nz_zip_factor, ncoils);
    end
    
    % X dimension, fft & remove oversampling
    kx1 = nx_recon2/nx_oversample;
    sub_datax = fftshift(ifft(sub_data,nx_recon2*nx_oversample,1),1);
    sub_datax = sub_datax(kx1+1:kx1+nx_recon2,:,:,:);
    
    ctrl_datax = fftshift(ifft(ctrl_data,nx_recon2*nx_oversample,1),1);
    ctrl_datax = ctrl_datax(kx1+1:kx1+nx_recon2,:,:,:);
    
    if (use_M0)
        M0_datax = fftshift(ifft(M0_data,nx_recon2*nx_oversample,1),1);
        M0_datax = M0_datax(kx1+1:kx1+nx_recon2,:,:);
    end
    
    if (debug_mode)
        lbl_datax = fftshift(ifft(lbl_data,nx_recon2*nx_oversample,1),1);
        lbl_datax = lbl_datax(kx1+1:kx1+nx_recon2,:,:,:);
    end
    
    % Peform Z fft, if 3D
    if (nz > 1)
        sub_datax = ifft(sub_datax, nz*nz_zip_factor, 3);
        ctrl_datax = ifft(ctrl_datax, nz*nz_zip_factor, 3);
        if (debug_mode)
            lbl_datax = ifft(lbl_datax, nz*nz_zip_factor, 3);
        end
        if (use_M0)
            M0_datax = ifft(M0_datax, nz*nz_zip_factor, 3);
        end
    end
    
    % Rotate and perform homodyne along Y
    clear img;
    clear M0_img;
    hd_flag = 1;
    for (coil_ct = 1:ncoils)
        for (slice_ct = 1:nz*nz_zip_factor)
            %         [sub_img1, temp] = part_echo(sub_datax(:,:,slice_ct,coil_ct)', ny_recon, zip_factor*nx, zip_factor, hd_flag);
            [sub_img1, temp] = part_echo(sub_datax(:,:,slice_ct,coil_ct)', ny_recon, nx_recon2, zip_factor_ny, hd_flag); % Yue
            % img(:,:,slice_ct,coil_ct) = sub_img1'.*chopx.*chopy;
            img(:,:,slice_ct,coil_ct) = sub_img1';
            
            %         [ctrl_img1, temp] = part_echo(ctrl_datax(:,:,slice_ct,coil_ct)', ny_recon, zip_factor*nx, zip_factor, hd_flag);
            [ctrl_img1, temp] = part_echo(ctrl_datax(:,:,slice_ct,coil_ct)', ny_recon, nx_recon2, zip_factor_ny, hd_flag); % Yue
            % low_res_img(:,:,slice_ct,coil_ct) = temp'.*chopx.*chopy;
            low_res_img(:,:,slice_ct,coil_ct) = temp';
            
            if (debug_mode)
                ctrl_img(:,:,slice_ct,coil_ct) = ctrl_img1'.*chopx.*chopy;
%                             [lbl_img1, temp] = part_echo(lbl_datax(:,:,slice_ct,coil_ct)', ny_recon, zip_factor*nx, zip_factor, hd_flag);
                [lbl_img1, temp] = part_echo(lbl_datax(:,:,slice_ct,coil_ct)', ny_recon, nx_recon2, zip_factor_ny, hd_flag); % Yue
                % lbl_img(:,:,slice_ct,coil_ct) = lbl_img1'.*chopx.*chopy;
                lbl_img(:,:,slice_ct,coil_ct) = lbl_img1';%LZ
            end
        end
        
        % M0 image
        if (use_M0)
            %         [M0_img1, temp] = part_echo(M0_datax(:,:,coil_ct)', ny_recon, zip_factor*nx, zip_factor, hd_flag);
            [M0_img1, temp] = part_echo(M0_datax(:,:,coil_ct)', ny_recon_M0, nx_recon2, zip_factor_ny_M0, hd_flag); % Yue
            % M0_img(:,:,coil_ct) = M0_img1'.*chopx.*chopy;
            M0_img(:,:,coil_ct) = M0_img1';
            if (nz == 1)
                % low_res_img(:,:,1,coil_ct) = temp'.*chopx.*chopy;
                low_res_img(:,:,1,coil_ct) = temp';
            end
        end
        
    end
    
    fprintf('Reconstruction Done\n');
    
    % Create low-res coil sensitivity
    % coil_sens = zeros(nx_recon2, ny_recon2, nz*nz_zip_factor, ncoils); % Yue
    low_SS = sqrt(sum(abs(low_res_img).^2, 4));
    
    for (coil_ct = 1:ncoils)
        for (slice_ct = 1:nz*nz_zip_factor)
            coilSens = low_res_img(:,:,slice_ct,coil_ct);
            coil_sens(:,:,slice_ct,coil_ct) = (coilSens.*conj(coilSens)./(abs(coilSens) + 1e-10))./low_SS(:,:,slice_ct);
        end
    end
    
    % Create low-res coil sensitivity weighted coil combined image
    hd_flag1 = 1;
    if (hd_flag1)
        final_image = real( sum(img.*coil_sens,4));
        if (use_M0)
            M0_final_image = real( sum(M0_img.*squeeze(coil_sens),3));
        end
        
        if (debug_mode)
            final_ctrl_image = real( sum(ctrl_img.*coil_sens,4));
            final_lbl_image = real( sum(lbl_img.*coil_sens,4));
        end
    else
        %     final_image = sqrt( sum(img.*coil_sens,4));
        final_image = sqrt( sum(abs(img).*abs(img),4));
        if (use_M0)
            % M0_final_image = sqrt( sum(M0_img.*squeeze(coil_sens),3));
            M0_final_image = sqrt( sum(abs(M0_img).*abs(M0_img),3));
        end
        
        if (debug_mode)
            final_ctrl_image = sqrt( sum(ctrl_img.*coil_sens,4));
            final_lbl_image = sqrt( sum(lbl_img.*coil_sens,4));
        end
    end
    %take the absolute of the final image for Hd = 1 & Hd = 2
    if hd_flag1 == 1 || hd_flag1 == 2
        final_image = abs(final_image);
    end
    savefile = "../MATLAB/Images/ImageFiles/Raw005_2D_T2Co80_Hd1_Dec1_VarTest.mat";
    save(savefile,'final_image');
    figure();imshow(final_image,[])
    title("Final Averaged Subtraction Image");

    
else % Fully acquired kspace - only FFT needed
    
%      control = squeeze(data(:,:,1,:));
%     label = squeeze(data(:,:,2,:));% LZ
    
    control = squeeze(data(:,:,:,2,:));
    label = squeeze(data(:,:,:,1,:));
    sub = squeeze(control-label);
 
    if use_M0
        M0 = M0_data;
        
    end
    
    recon_resolution_factor = 2;
    
    siz = [recon_resolution_factor*size(data,1) recon_resolution_factor*size(data,2) recon_resolution_factor*size(data,3)];
    
    nx = size(sub,1);
    ny = size(sub,2);
    nz = size(sub,3);
    ncoil = size(sub,4);
    
    if use_M0
        M0_siz = [recon_resolution_factor*size(M0_data,1) recon_resolution_factor*size(M0_data,2) recon_resolution_factor*size(M0_data,3)];
        M0_nx = size(M0,1);
        M0_ny = size(M0,2);
        M0_nz = size(M0,3);
        M0_ncoil = size(M0,4);
        
    end
    %% Filter
    
    windowA = 16;
    windowB = 16;
    fermi_filter = fermi2d(nx,ny,nx/2,ny/2,windowA,windowB);
    if use_M0
        M0_fermi_filter = fermi2d(M0_nx,M0_ny,M0_nx/2,M0_ny/2,windowA,windowB);
    end
    for z = 1:nz
        for coil = 1:ncoil
            sub(:,:,z,coil) = sub(:,:,z,coil).*fermi_filter;
            label(:,:,z,coil) = label(:,:,z,coil).*fermi_filter;
            control(:,:,z,coil) = control(:,:,z,coil).*fermi_filter;
%         sub(:,:,coil) = sub(:,:,coil).*fermi_filter;
%             label(:,:,coil) = label(:,:,coil).*fermi_filter;
%             control(:,:,coil) = control(:,:,coil).*fermi_filter;% LZ 
%             sub(:,:,z,coil) = sub(:,:,z,coil);
%             label(:,:,z,coil) = label(:,:,z,coil);
%             control(:,:,z,coil) = control(:,:,z,coil); %% without fermi filter
        end
    end
    
    
    
    for coil = 1:ncoil
        
        SUB_IM(:,:,:,coil) = ifftshift(ifftn(sub(:,:,:,coil),siz),1);
        CONTROL_IM(:,:,:,coil) = ifftshift(ifftn(control(:,:,:,coil),siz),1);
        LABEL_IM(:,:,:,coil) = ifftshift(ifftn(label(:,:,:,coil),siz),1);
        
%           SUB_IM(:,:,coil) = ifftshift(ifftn(sub(:,:,coil),siz),1);
%         CONTROL_IM(:,:,coil) = ifftshift(ifftn(control(:,:,coil),siz),1);
%         LABEL_IM(:,:,coil) = ifftshift(ifftn(label(:,:,coil),siz),1);% LZ
    end
    final_image = sqrt(sum(abs(squeeze(SUB_IM)).*abs(squeeze(SUB_IM)),4)); % [x y ndyn] double
    final_image = final_image(size(final_image,1)*.25:size(final_image,1)*.75,:,:); % Crop oversampled parts of image
    
    
    
    if use_M0
        
        for z = 1:nz
            for coil = 1:ncoil
                M0(:,:,z,coil) = M0(:,:,z,coil).*M0_fermi_filter;
%                 M0(:,:,z,coil) = M0(:,:,z,coil);
            end
        end
        
        
        for coil = 1:ncoil
            
            M0_IM(:,:,:,coil) = ifftshift(ifftn(M0(:,:,:,coil),M0_siz),1);
            
        end
        
        M0_final_image = sqrt(sum(abs(squeeze(M0_IM)).*abs(squeeze(M0_IM)),4)); % [x y ndyn] double
        M0_final_image = M0_final_image(size(M0_final_image,1)*.25:size(M0_final_image,1)*.75,:,:); % Crop oversampled parts of image
        
        
    end 
end

if nz>1
    implay(final_image)
    implay(M0_final_image)
end

%% Get M0 Mean value
% figure; imagesc(M0_final_image); colormap gray; title('Draw ROI + double-click')
% % figure; imagesc(M0_final_image(:,:,21)); colormap gray; title('Draw ROI + double-click')
% h=imellipse;
% position = wait(h); % double click ellipse when finished positioning
% BW = createMask(h);
% M0_roi=BW.*M0_final_image;
% M0_roi(M0_roi==0) =nan;
% M0_mean = nanmean(nanmean(M0_roi));

if (use_M0)
    lambda  = 0.9;      % Partition coefficient
    alpha   = 0.6;      % Labeling efficiency
    tau     = 1.5;      % Labeling time
    t_meas  = 3.0;      % label (1.5) + post-label delay (1.5)
    delta_t = 0.75;     % Arrival time
    T1blood = 1.6;      % 1.3 s for 1.5T and 1.6 s for 3T
    
    s1 = exp(-delta_t/T1blood);
    s2 = 1 - exp(-tau/T1blood);
    s3 = exp(-(t_meas-tau-delta_t)/T1blood);

    thresh1 = 0.05*max(final_image(:));
    perf_mask = (final_image>thresh1);

    max_M0  = max(M0_final_image(:));
    thresh2 = 0.1*max_M0;
    M0_mask = (M0_final_image>thresh2);
 
    mask_image = perf_mask.*thresh2;
% Detect different sized pCASL/M0 data
%M0_mask= 1;
% M0_mask = 1;
% mask_image =1;

if size(final_image) ~= size(M0_final_image)
    figure; imshow(M0_final_image,[]); title('M0 and ASL do not match. Draw Kidney ROI');
%     figure; imshow(M0_final_image(:,:,21),[]); title('M0 and ASL do not match. Draw Kidney ROI');
    mask = roipoly;
    M0_single_value = mean(M0_final_image(mask == 1));
    
    perflow_data = final_image./(M0_single_value); %if not same size, need to draw M0 ROI and divide by mean M0
    
else
    perflow_data = final_image./(M0_final_image); %if not same size, need to draw M0 ROI and divide by mean M0
end


    perflowdata = (perflow_data*lambda)/(2*alpha*T1blood*s1*s2*s3);
    perflow = M0_mask.*perflowdata*100*60;       % Scale to get units: ml/100ml.min
   
%     perflow_data = final_image./(M0_mean); %if not same size, need to draw M0 ROI and divide by mean M0
%     perflowdata = (perflow_data*lambda)/(2*alpha*T1blood*s1*s2*s3);
%     perflow_mean = M0_mask.*perflowdata*100*60;       % Scale to get units: ml/100ml.min

    figure; imshow(perflow,[]); colormap (gca,jet); colorbar;
end

for (slice_ct = 1:nz*nz_zip_factor)
    if (debug_mode)
        figure();
        subplot(221); imshow(final_ctrl_image(:,:,slice_ct),[]); title('control');
        subplot(222); imshow(final_lbl_image(:,:,slice_ct),[]); title('label');
        subplot(223); imshow(final_image(:,:,slice_ct),[]); title('pcasl');
        if (use_M0)
            subplot(224); imshow(perflow(:,:,slice_ct),[0 600]); title('flow');
        end
%         title(sprintf('slice #%d of %d', slice_ct, nz));
        
    else
        figure();
        
        if (use_M0)
            %             subplot(221); imshow(final_image(:,:,slice_ct),[]); title('\Delta M')
            %             subplot(222); imshow(M0_final_image(:,:,slice_ct),[]); title('M_0')
            %             subplot(223); imshow(perflow(:,:,slice_ct),[0 600]); title('Perfusion')
            figure; imshow(final_image(:,:,slice_ct),[]); title('\Delta M')
            figure; imshow(M0_final_image(:,:,slice_ct),[]); title('M_0')
            figure; imshow(perflow(:,:,slice_ct),[0 600]); title('Perfusion')
        else
            imshow(final_image(:,:,slice_ct),[]);
        end
        title(sprintf('slice #%d of %d', slice_ct, nz*nz_zip_factor));
        if (nz > 1)
            pause;
        end
    end
end

if (write_data)
    
    %     final_image = sum(img, 4);
    % uigetdir()
    [fname_dicomread, dir_name_dicom] = uigetfile('*.*', 'Choose dicom file',dir_name); %DD adds
   % uigetfile('C:\Users\S172551\Downloads\MRI data\');
   % dir_name_dicom = 'C:\Users\S172551\Downloads\MRI data\20161104_RCC_093\DICOM';
    
  %  fname_dicomread = 'IM_0025';
    
    min_val = min(min(final_image(:)));
    max_val = max(max(final_image(:)));
    
    for (slice_ct = 1:nz*nz_zip_factor)
        
        diff_dicom_image = floor(4096*(final_image(:,:,slice_ct) - min_val)/(max_val - min_val));
        % Yue: 0 filling in image space
        %         diff_dicom_image_full = zeros(ny_recon2,ny_recon2);
        %         diff_dicom_image_full((ny_recon2-nx_recon2)/2+1:(ny_recon2-(ny_recon2-nx_recon2)/2), :) = diff_dicom_image;
        diff_dicom_image_full = padzerosoutside2d(diff_dicom_image);
        
        % Read DICOM data
        dicom_filename = [dir_name_dicom '\' fname_dicomread];
        dicom_info = dicominfo(dicom_filename);
        dicom_data = dicomread(dicom_filename);
        
        dicom_info.SeriesNumber = 100*series_no+12; % Yue: 12 for List Data recon
        dicom_info.SeriesDescription = '2D Complex Difference from Data List';
        dicom_info.ProtocolName = '2D PCASL Complex Difference from Data List';
        if (nz > 1) % 3D
            dicom_info.SeriesDescription = '3DFSE Complex Difference';
            dicom_info.ProtocolName = '3D PCASL Complex Difference';
        end
        %         dicom_info = rmfield(dicom_info,'PrivateInformationCreatorUID'); %Yue: Pay attention that not all dicom has PrivateInformationCreatorUID
        
        %         dicom_info.SeriesNumber = 200*series_no+11;
        %         dicom_info.SeriesDescription = 'PCASL Sum of Coils';
        %         dicom_info.ProtocolName = 'PCASL Sum of Coils';
        
        dicom_info.InstanceNumber = slice_ct;
        if (nz > 1)
            dicom_info.Private_2005_1008 = dicom_info.Private_2005_1008 + (slice_ct-1)*slice_thickness;
            if (sagittal)
                dicom_info.ImagePositionPatient(1) = dicom_info.ImagePositionPatient(1) + (slice_ct-1)*slice_thickness;
            else
                dicom_info.ImagePositionPatient(2) = dicom_info.ImagePositionPatient(2) + (slice_ct-1)*slice_thickness;
            end
        end
        
        dicom_info.WindowCenter = (4096 + 0)/2;
        dicom_info.WindowWidth = (4096 - 0);
        
        %         dicom_info.Rows = nx*zip_factor;
        %         dicom_info.Columns = ny_recon*zip_factor;
        %         dicom_info.Rows = nx_recon2; % Yue
        %         dicom_info.Columns = ny_recon2; %Yue
        dicom_info.Width = dicom_info.Columns;
        dicom_info.Height = dicom_info.Rows;
        if (nz > 1)
            row_spacing = dicom_info.ReconstructionDiameter/dicom_info.Rows;
            col_spacing = dicom_info.ReconstructionDiameter/dicom_info.Columns;
            dicom_info.PixelSpacing = [row_spacing; col_spacing];
        end
        
        diff_dicomname = sprintf('se%d_fromListData_img%d.dcm',dicom_info.SeriesNumber,slice_ct);% Yue
        diff_dicom_filename = [dir_name_dicom diff_dicomname];
        dicomwrite(uint16(diff_dicom_image_full), diff_dicom_filename, dicom_info, 'CreateMode', 'copy');
    end
end

% Write Perfusion Quantified data
if (write_quant_data)
    min_val = min(min(perflow(:)));
    max_val = max(max(perflow(:)));
    
    for (slice_ct = 1:nz*nz_zip_factor)
        
        perflow_dicom_image = floor(perflow(:,:,slice_ct));
        % Yue: 0 filling in image space
        %         perflow_dicom_image_full = zeros(ny_recon2,ny_recon2);
        %         perflow_dicom_image_full((ny_recon2-nx_recon2)/2+1:(ny_recon2-(ny_recon2-nx_recon2)/2), :) = perflow_dicom_image;
        perflow_dicom_image_full = padzerosoutside2d(perflow_dicom_image);
        save('perflow.mat','perflow_dicom_image_full');
        % Read DICOM data
        dicom_filename = [dir_name_dicom  '\' fname_dicomread];
        dicom_info = dicominfo(dicom_filename);
        dicom_data = dicomread(dicom_filename);
        
        dicom_info.SeriesNumber = 100*series_no+13; % Yue: 13 for List Data Perfusion
        dicom_info.SeriesDescription = '2D Perfusion Map from Data List after Scaling';
        dicom_info.ProtocolName = '2D PCASL Perfusion Map from Data List after Scaling';
        dicom_info.InstanceNumber = slice_ct;
        
        
        % check if we need to adjust rescaling
        % imshow
        if isfield(dicom_info,'PerFrameFunctionalGroupsSequence')
            
            % Yue: Correct the scalling issue
%             for scale_count = 1 : 2 * ndynamics
            for scale_count = 1
                item_name = sprintf('Item_%d', scale_count);
                if dicom_info.PerFrameFunctionalGroupsSequence.(item_name).PixelValueTransformationSequence.Item_1.RescaleSlope ~= 1
                    dicom_info.PerFrameFunctionalGroupsSequence.(item_name).PixelValueTransformationSequence.Item_1.RescaleSlope = 1.0;
                end
                if dicom_info.PerFrameFunctionalGroupsSequence.(item_name).PixelValueTransformationSequence.Item_1.RescaleIntercept ~= 0
                    dicom_info.PerFrameFunctionalGroupsSequence.(item_name).PixelValueTransformationSequence.Item_1.RescaleIntercept = 0;
                end
                % dicom_info.PerFrameFunctionalGroupsSequence.Item_1.PixelValueTransformationSequence.Item_1.RescaleIntercept
                % RI: 0028, 1052 (Resale Intercept) Original value = 0;
                % RS: 0028, 1053 (Rescale Slope) Original value = 3.626862026862020;
                %         end
            end
        end
        
        % Yue: Correct the scalling issue end
        if (nz > 1)
            dicom_info.Private_2005_1008 = dicom_info.Private_2005_1008 + (slice_ct-1)*slice_thickness;
            dicom_info.ImagePositionPatient(2) = dicom_info.ImagePositionPatient(2) + (slice_ct-1)*slice_thickness;
        end
        dicom_info.WindowCenter = (min_val + max_val)/2;
        dicom_info.WindowWidth = (max_val - min_val);
        
        %         dicom_info.Rows = nx*zip_factor;
        %         dicom_info.Columns = ny_recon*zip_factor;
        %         dicom_info.Rows = nx_recon2; % Yue
        %         dicom_info.Columns = ny_recon2; %Yue
        dicom_info.Width = dicom_info.Columns;
        dicom_info.Height = dicom_info.Rows;
        
        %         dicom_info = rmfield(dicom_info,'PrivateInformationCreatorUID'); %Yue: Pay attention that not all dicom has PrivateInformationCreatorUID
        %         row_spacing = dicom_info.ReconstructionDiameter/dicom_info.Rows;
        %         col_spacing = dicom_info.ReconstructionDiameter/dicom_info.Columns;
        %         dicom_info.PixelSpacing = [row_spacing; col_spacing];
        
        quantfname = sprintf('se%d_fromListData_img%d_quant_after_scaling.dcm',dicom_info.SeriesNumber,slice_ct); % Yue
        newquantdicom_filename = [dir_name_dicom quantfname];
        dicomwrite(uint16(perflow_dicom_image_full), newquantdicom_filename, dicom_info, 'CreateMode', 'copy');
    end
end
 fprintf('Export Dicom\n');
 toc