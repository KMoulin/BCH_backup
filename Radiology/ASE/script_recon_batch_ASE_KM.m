function script_recon_batch_ASE_KM(struct_diff)
table=[];
Dcm=[];
listing=struct_diff.Listing;
for cpt=1:1:length(listing)
        [FolderName, name, fExt] = fileparts(fullfile(listing(cpt).folder , listing(cpt).name));
        
        if (strcmp(fExt, '.dcm') | strcmp(fExt, '.IMA') | isempty(fExt)) & ~listing(cpt).isdir % if listing(cpt).name(end-2:end) == 'dcm' | listing(cpt).name(end-2:end) == 'IMA'
            
                Dcm(:,:,:,cpt)=double(dicomread(fullfile(FolderName , listing(cpt).name)));
                info=dicm_hdr(fullfile(FolderName , listing(cpt).name));               
                [table(:,cpt), names,units]=comment_extract_local(info.ImageComments);
                               
                if isfield(info,'PerFrameFunctionalGroupsSequence')
                    % Temporal_pos(cpt)=info.PerFrameFunctionalGroupsSequence.Item_1.FrameContentSequence.Item_1.TemporalPositionIndex;
                    %Echo_n(cpt)=info.PerFrameFunctionalGroupsSequence.Item_1.MREchoSequence.Item_1.EchoNumber;
                    Echo(cpt)=info.PerFrameFunctionalGroupsSequence.Item_1.MREchoSequence.Item_1.EffectiveEchoTime;
                else
                    Echo(cpt)=info.EchoTime;     
                end
        end
    end

    [S0map, T2map, R2primeMap, fitInfo] = fitASE_R2prime(Dcm, Echo, table(1,:)*1e-3);
    save(fullfile(struct_diff.ReconFolder ,'T2map.mat'),'S0map','T2map','R2primeMap','fitInfo');
end
function [values,names,units]=comment_extract_local(str)
% Example DICOM image comment string
%str = "Not for diagnostic use, e2, DTE = 3000 us, ZMoment = 0e-6 T*ms/m, Rep = 75, Echo = 1, MaxDTE = 44000 us, dTE = 1000us";

% Split on commas
parts = strtrim(split(str, ','));

% Prepare containers
names  = {};
values = [];
units  = {};

for i = 1:numel(parts)
    p = parts{i};
    
    % Match pattern: name = value unit  (unit optional, spacing optional)
    tok = regexp(p, '^(\w+)\s*=\s*([-+]?[\d.eE+-]+)\s*(\S*)$', 'tokens');
    
    if ~isempty(tok)
        names{end+1}  = tok{1}{1};                 %#ok<*AGROW>
        values(end+1) = str2double(tok{1}{2});
        units{end+1}  = tok{1}{3};
    end
end

% Build table
% paramTable = table(names', values', units', ...
%     'VariableNames', {'Parameter','Value','Unit'});
% 
% disp(paramTable)
end
function [S0map, T2map, R2primeMap, fitInfo] = fitASE_R2prime(imgStack, TE_ms, DTE_ms, varargin)
% FITASE_R2PRIME  Voxelwise linear fit of ASE data to extract S0, T2, R2'
%
%   [S0map, T2map, R2primeMap, fitInfo] = fitASE_R2prime(imgStack, TE_ms, DTE_ms)
%
%   INPUTS
%     imgStack : [rows x cols x (slices) x N] image stack, N = number of
%                acquisitions (all echoes x all DTE steps stacked together)
%     TE_ms    : [N x 1] echo time (ms) for each acquisition
%     DTE_ms   : [N x 1] DTE offset (ms, signed) for each acquisition
%
%   NAME-VALUE OPTIONS
%     'tc_ms'  : exclusion half-width around DTE=0 (default 1.5 ms)
%     'mask'   : logical mask [rows x cols x (slices)] restricting the fit
%
%   OUTPUTS
%     S0map, T2map, R2primeMap : parameter maps (same spatial size as imgStack)
%     fitInfo                  : struct with residuals, R^2, included points

    p = inputParser;
    addParameter(p, 'tc_ms', 1.5);
    addParameter(p, 'mask', []);
    parse(p, varargin{:});
    tc_ms = p.Results.tc_ms;
    mask  = p.Results.mask;

    sz = size(imgStack);
    N  = sz(end);
    spatialSize = sz(1:end-1);

    TE_ms  = TE_ms(:);
    DTE_ms = DTE_ms(:);
    assert(numel(TE_ms) == N && numel(DTE_ms) == N, ...
        'TE_ms and DTE_ms must have one entry per image in imgStack');

    % --- exclude near-zero DTE points (non-linear "rounding" regime) ---
    keep = abs(DTE_ms) >= tc_ms;
    fprintf('Excluding %d/%d points with |DTE| < %.2f ms\n', ...
        sum(~keep), N, tc_ms);

    TEk  = TE_ms(keep);
    DTEk = DTE_ms(keep);
    Nk   = numel(TEk);

    % --- reshape image stack to [N x Nvoxels] ---
    imgReshaped = reshape(imgStack, [], N)';     % N x Nvoxels
    imgReshaped = imgReshaped(keep, :);          % Nk x Nvoxels

    if ~isempty(mask)
        voxelIdx = find(mask(:));
    else
        voxelIdx = find(all(imgReshaped > 0, 1)); % avoid log(0) / negative noise
    end

    Y = imgReshaped(:, voxelIdx);
    Y(Y <= 0) = eps;                             % guard against non-positive signal
    logY = log(Y);                               % Nk x Nvoxels_valid

    % --- design matrix: ln S = a - TE*b - |DTE|*c ---
    A = [ones(Nk,1), -TEk, -abs(DTEk)];           % Nk x 3

    % --- solve all voxels simultaneously ---
    coeffs = A \ logY;                            % 3 x Nvoxels_valid

    a = coeffs(1,:);
    b = coeffs(2,:);
    c = coeffs(3,:);

    S0_valid       = exp(a);
    T2_valid       = 1 ./ b;
    R2prime_valid  = c;

    % --- residuals / goodness of fit ---
    logY_hat = A * coeffs;
    resid    = logY - logY_hat;
    SSres    = sum(resid.^2, 1);
    SStot    = sum((logY - mean(logY,1)).^2, 1);
    R2_valid = 1 - SSres ./ SStot;

    % --- scatter results back into full image grid ---
    S0map      = zeros(prod(spatialSize),1);
    T2map      = zeros(prod(spatialSize),1);
    R2primeMap = zeros(prod(spatialSize),1);
    R2fitMap   = zeros(prod(spatialSize),1);

    S0map(voxelIdx)      = S0_valid;
    T2map(voxelIdx)       = T2_valid;
    R2primeMap(voxelIdx) = R2prime_valid;
    R2fitMap(voxelIdx)   = R2_valid;

    S0map      = reshape(S0map, spatialSize);
    T2map      = reshape(T2map, spatialSize);
    R2primeMap = reshape(R2primeMap, spatialSize);
    R2fitMap   = reshape(R2fitMap, spatialSize);

    fitInfo.R2fitMap   = R2fitMap;
    fitInfo.nIncluded  = Nk;
    fitInfo.nExcluded  = N - Nk;
    fitInfo.tc_ms      = tc_ms;
end
