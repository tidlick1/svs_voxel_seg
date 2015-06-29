function svs_voxel_seg(root_dir)

%Developed by Scott Quadrelli June 2015
% The aim of this function is to perform automated voxel partial volume
% segmentation for a collection of .rda files
% This program makes use of GannetMask_Siemens function, please reference
% accordingly (Harris, AD, Puts, NAJ, Edden, RAE. 2015.)

% folder = the directory where the folder of structural image and .rda
% files are contained.

global cmdout;
global cmdout_nonpar;
fileID = fopen([root_dir '/results.txt'],'w');
addpath /Users/quadrellis/spm12;
setenv('Data',root_dir);

%check FSL installation
global fsldir;
fsldir = '/usr/local/fsl/';
setenv('FSLDIR',fsldir);  % this to tell where FSL folder is
if ~exist(fsldir,'dir')
	error('%s: error fsldir (%s) not found',mfilename, fsldir);
end
setenv('PATH',sprintf('%s:%s',fullfile(fsldir,'bin'),getenv('PATH')));
bet = [fsldir '/bin/bet'];
if ~exist(bet)
	error('%s: error %s not found',mfilename,bet);
end

%determines the number of directories in root
subdirs = dir(root_dir);
%This filters out all the items in the main folder that are not directories
subdirs(~[subdirs.isdir]) = [];
%And this filters out the parent and current directory '.' and '..'
tf = ismember( {subdirs.name}, {'.', '..'});
subdirs(tf) = [];
no_folders = length(subdirs);

for i = 1: no_folders
    currentdir = subdirs(i).name;
    folder = [root_dir '/' currentdir];
    %finds the rda and the structural
    rdafiles=dir(sprintf('%s/*.rda',folder)); 
    structuralfile=dir(sprintf('%s/*.nii',folder));
    %Creates folders for the outputs
    if ~exist([folder '/structural'],'dir')
        mkdir(sprintf('%s/structural',folder)); %brain extraction output
    end
    if ~exist((sprintf('%s/fast_output',folder)),'dir')
        mkdir(sprintf('%s/fast_output',folder)); %segmentation output
    end
    if ~exist((sprintf('%s/masks',folder)),'dir')
        mkdir(sprintf('%s/masks',folder)); %voxel_masks
    end
    
    structuralout = (sprintf('%s/structural',folder));
    segout = (sprintf('%s/fast_output',folder));
    structural_location = (sprintf('%s/%s',folder, structuralfile.name));

    %performs brain extraction, after confiming it hasn't already been done
    if ~exist([structuralout '/extractedbrain.nii.gz'],'file')
        disp(['Running brain extraction on' structuralfile.name])
        fslCmd(sprintf('bet %s/%s %s/extractedbrain -f 0.5 -g 0', folder,(structuralfile.name), structuralout));
    else
        disp(['Brain extraction has already been done for' structuralfile.name])
    end
    
    %Partial Volume Segmentation, , after confiming it hasn't already been done
    if ~exist([segout '/fast_output_seg.nii.gz'],'file')
        disp(['Running brain segmentation on' structuralfile.name])
        fslCmd(sprintf('fast -t 1 -n 3 -H 0.1 -I 4 -l 20.0 -o %s/fast_output %s/extractedbrain', segout, structuralout));
    else
        disp(['Brain segmentation has already been done for' structuralfile.name])
    end
    
    %Creates voxel masks for each .rda and saves in the masks folder
        for file = rdafiles'
            currentrda = (sprintf('%s/%s',folder,file.name));
            GannetMask_Siemens(currentrda,structural_location);
              %call fsl stats to determine fractions
              %Uses regular expressions to just get the end of the mask filename
              currentmask = regexp(currentrda, '/[\t 0-9A-Z_a-z]+(\.rda)$','match');
              %Strips off the .rda section of the filename
              currentmask = strrep(currentmask,'.rda','');
              fslCmd_nonpar(sprintf('fslstats -t %s/fast_output_pve_0 -k %s/masks%s_mask.nii -m',segout,folder,char(currentmask)));
              greyfrac = cmdout_nonpar;
              fslCmd_nonpar(sprintf('fslstats -t %s/fast_output_pve_1 -k %s/masks%s_mask.nii -m',segout,folder,char(currentmask)));
              whitefrac = cmdout_nonpar;
              fslCmd_nonpar(sprintf('fslstats -t %s/fast_output_pve_2 -k %s/masks%s_mask.nii -m',segout,folder,char(currentmask)));
              csffrac = cmdout_nonpar;
              %Combines each fraction into an array
              tot_fractions= [currentdir,currentmask,cellstr(greyfrac),cellstr(whitefrac),cellstr(csffrac)]
              %outputs the fractions to a textfile called results, comma seperated
              fprintf(fileID, '%s,%s,%s,%s,%s\n', tot_fractions{:});      
        end 
end

fclose(fileID);
clear all


function fslCmd (Cmd);
%execute a fsl command, e.g. fslCmd('fslinfo a.nii');
global fsldir;
global cmdout;
        %spmd
        command=sprintf('sh -c ". %setc/fslconf/fsl.sh; %sbin/%s"\n',fsldir,fsldir, Cmd);
        [status,cmdout] = system(command);
    %end
%end fslCmd()

function fslCmd_nonpar (Cmd);
%execute a fsl command not in parallel;
global fsldir;
global cmdout_nonpar;
        command=sprintf('sh -c ". %setc/fslconf/fsl.sh; %sbin/%s"\n',fsldir,fsldir, Cmd);
        [status_nonpar,cmdout_nonpar] = system(command);