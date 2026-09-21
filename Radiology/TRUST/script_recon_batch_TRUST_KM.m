function script_recon_batch_TRUST_KM(struct_diff,GUI)
    
    % Batch script meant for being called from the Batch_Manager. 
    % the script expect a struct_diff 
    warning off;

    struct_diff.ReconFolder=char(struct_diff.ReconFolder);
    mkdir(struct_diff.ReconFolder);

    enum=[];
    Dcm=[];
    
  
    enum.dcm_dir=char(struct_diff.DcmFolder);
    enum.recon_dir=char(struct_diff.ReconFolder);

    [Dcm,enum2]= DiffRecon_ToolBox.Read_all_KM(struct_diff.Listing); %Read_all_4DFlow_KM

    csa = dicm_hdr(struct_diff.Listing(1).name);
    save(fullfile(enum.recon_dir ,'RAWn.mat'),'Dcm','enum2','csa');

    %%
    tmp_trust_acq=string([]);
    tmp_trust_avg=[];
    
    tmp_num=string([]);
    
    tmp_acq=[];
    tmp_avg=[];
    tmp_time=[];
    Dcm2=[];
    
    % Step 1: We organize the dicoms by acq and avg
    for cpt=1:1:length(enum2.header)
    
        tmp_name=enum2.header(cpt).dicomheader.ProtocolName;
        %disp(tmp_name)
        tmp_num(cpt)=tmp_name(end);
        
        if isempty(tmp_trust_acq)
            tmp_trust_acq(1)=tmp_num(cpt);
            tmp_trust_avg(1)=1;
        else
            if ~contains(tmp_trust_acq,tmp_num(cpt)) 
                tmp_trust_acq(end+1)=tmp_num(cpt);
                tmp_trust_avg(end+1)=1;
            else
                idx=find(strcmp(tmp_trust_acq,tmp_num(cpt)));
                tmp_trust_avg(idx)=tmp_trust_avg(idx)+1;    
            end
        end
        idx=find(strcmp(tmp_trust_acq,tmp_num(cpt)));
        tmp_acq(cpt)=idx;
        tmp_avg(cpt)=tmp_trust_avg(idx);
        tmp_time(tmp_acq(cpt),tmp_avg(cpt))=str2num(enum2.header(cpt).dicomheader.ContentTime);
        Dcm2(:,:,1,tmp_acq(cpt),tmp_avg(cpt))=Dcm(:,:,cpt);
    end
   
    % Step2:  We sort the acquisitions by scan time
    
    for cpt=1:1:size(tmp_time,1)
        tmp_time_sort=tmp_time(cpt,:);
        [~,idx]=sort(tmp_time_sort);
        Dcm2(:,:,1,cpt,:)=Dcm2(:,:,1,cpt,idx);
    end
%     save(fullfile(enum.recon_dir ,'RAW_TRUST.mat'),'Dcm2','enum2');
    
    % Step3: Calculate the difference from SAT reference
    Dcm3=mean(Dcm2(:,:,:,:,2:2:end),5)-mean(Dcm2(:,:,:,:,1:2:end),5);
    mip_mask=max(Dcm3,[],4);
    mip_mask(mip_mask<20)=nan;
    mip_mask(~isnan(mip_mask))=1;
    
%    save(fullfile(enum.recon_dir ,'Mean_TRUST.mat'),'Dcm3','mip_mask','enum2');

    % Step4: T2fitting
    listTEs=str2double(tmp_trust_acq)*20; %ms;
    [T2map, S0map]=fitT2_local(Dcm3,listTEs,mip_mask);
    
    %Oxmap=T22Ox_local(T2map);
    Oxmap = estimate_oxygenation(T2map*1e-3);
    save(fullfile(enum.recon_dir ,'T2_OX_TRUST.mat'),'T2map','S0map','Oxmap','enum2','listTEs');

function [T2map, S0map]=fitT2_local(Dcm,eTEs,mask)
% --- Inputs --------------------------------------------------------------
% diff_images: 4D array [x, y, nSlices, nETE]
% eTEs: vector of effective echo times in ms, e.g. [0 20 40 60 80]

[nx, ny, nz, nETE] = size(Dcm);

% --- Output maps ---------------------------------------------------------
T2map = zeros(nx, ny, nz);
S0map = zeros(nx, ny, nz);

% --- Fit options ---------------------------------------------------------
opts = optimset('Display', 'off');

% Monoexponential model: S = S0 * exp(-eTE / T2)
model = @(p, eTE) p(1) * exp(-eTE / p(2));

% Initial guess and bounds [S0, T2]
p0     = [1000, 80];   % initial guess
lb     = [0,     1];   % lower bounds
ub     = [1e6, 500];   % upper bounds (cap T2 at 500 ms to avoid CSF blowup)

% --- Voxelwise fitting ---------------------------------------------------
for iz = 1:nz
    for iy = 1:ny
        for ix = 1:nx
            if mask(ix,iy,iz)==1
            signal = squeeze(Dcm(ix, iy, iz, :));  % [nETE x 1]
         
            % Fit
            try
                p = lsqcurvefit(model, p0, eTEs(:), signal(:), lb, ub, opts);
                T2map(ix, iy, iz) = p(2);
                S0map(ix, iy, iz) = p(1);
            catch
                % Leave as zero if fit fails
            end
            end

        end
    end
end

S0map(S0map==ub(1))=0;
T2map(T2map==ub(2))=0;


end

function Oxmap=T22Ox_local(T2map)
    % Constants at 3T
A = 0.000601;   % ms^-1
B = 0.000436;   % ms^-1
C = 0.000557;   % ms^-1

% Assumed or measured hematocrit (typical assumed values)
Hct = 0.44;     % male:   ~0.44
% Hct = 0.40;   % female: ~0.40

% --- From a single T2 value (e.g. from ROI in sinus) --------------------
T2_blood = 85;  % ms, from your monoexponential fit

R2 = 1 / T2_blood;   % convert to R2 (ms^-1)

% Solve for Yv
% R2 = A + B*Hct + C*Hct*(1-Yv)^2
% (1-Yv)^2 = (R2 - A - B*Hct) / (C*Hct)
% Yv = 1 - sqrt(...)

Yv = 1 - sqrt((R2 - A - B*Hct) / (C * Hct));

% --- From a T2 map (voxelwise) -------------------------------------------
% T2map: [nx, ny, nz] in ms

R2map = 1 ./ T2map;                          % convert to R2
R2map(T2map == 0) = 0;                        % handle masked voxels

inner  = (R2map - A - B*Hct) / (C * Hct);
inner  = max(inner, 0);                       % clamp negative values (unphysical)

Oxmap = 1 - sqrt(inner);


end
function Yn = estimate_oxygenation(T2)
% ESTIMATE_OXYGENATION  Extracts oxygen saturation from blood T2
% using the Lu et al. MRM 2012 model:
%   1/T2 = A + B*(1-Yn) + C*(1-Yn)^2
%
% INPUTS:
%   T2  - Measured blood T2 in seconds (scalar or array)
%   A   - Model coefficient (s^-1)
%   B   - Model coefficient (s^-1)
%   C   - Model coefficient (s^-1)
%
% OUTPUT:
%   Yn  - Oxygen saturation (0 to 1)
%
% Example:
%   Yn = estimate_oxygenation(0.080, A, B, C)
    % Value in s^-1 for a CPMG of 5ms from Lu et al. MRM 2012  10.1002/mrm.22970 
    Hct      = 0.44;
    a1=-4.4;
    a2=39.1;
    a3=-33.5;
    b1=1.5;
    b2=4.7;
    c1=167.8;
    A =a1+ a2*Hct + a3 *Hct*Hct; 
    B= b1*Hct + b2*Hct*Hct;
    C= c1*Hct*(1-Hct);

    R2 = 1 ./ T2;  % Convert to R2 (s^-1)
    
    % Quadratic: C*x^2 + B*x + (A - R2) = 0
    % where x = (1 - Yn)
    discriminant = B^2 - 4*C*(A - R2);
    
    % Check for invalid values
    if any(discriminant(:) < 0)
        warning('Negative discriminant detected — T2 may be outside valid model range.');
    end
    
    % Take the physically meaningful root (smaller x)
    x = (-B + sqrt(discriminant)) / (2*C);
    
    % Recover Yn
    Yn = 1 - x;
    
    % Clip to physiological range [0, 1] and warn
    if any(Yn(:) < 0 | Yn(:) > 1)
        warning('Yn out of [0,1] range — check T2 values and coefficients.');
    end
    Yn = max(0, min(1, Yn));

end


end



