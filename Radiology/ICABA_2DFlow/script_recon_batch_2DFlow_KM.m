function script_recon_batch_2DFlow_KM(struct_diff,GUI)
    
    % Batch script meant for being called from the Batch_Manager. 
 
    warning off;

    struct_diff.ReconFolder=char(struct_diff.ReconFolder);
    mkdir(struct_diff.ReconFolder);

    enum=[];
    Dcm=[];
    
  
    enum.dcm_dir=char(struct_diff.DcmFolder);
    enum.recon_dir=char(struct_diff.ReconFolder);

    [Dcm,enum2]= DiffRecon_ToolBox.Read_all_KM(struct_diff.Listing); %Read_all_4DFlow_KM

    csa = dicm_hdr(struct_diff.Listing(1).name);

    tmp_name=enum2.header(1).dicomheader.ImageType;
    tmp_type=split(tmp_name,'\');
        if strcmp(tmp_type{3},'M')
            save(fullfile(enum.recon_dir ,'RAWn_M.mat'),'Dcm','enum2','csa');           
        elseif strcmp(tmp_type{3},'P')

            tmp_recon=split(struct_diff.ReconFolder,'_');
            
            tmp_recon{end-1}=num2str((str2num(tmp_recon{end-1})-2));
            tmp_recon{end-2}='';
            for cpt_str=1:1:length(tmp_recon)
                if cpt_str==1
                    path2mask=tmp_recon{cpt_str};
                else
                    if ~isempty(tmp_recon{cpt_str})
                        path2mask=append(path2mask, "_", tmp_recon{cpt_str});
                    end
                end
            end
            save(fullfile(enum.recon_dir ,'RAWn_P.mat'),'Dcm','enum2','csa','path2mask');
        else
            save(fullfile(enum.recon_dir ,'RAWn.mat'),'Dcm','enum2','csa');
        end
    end