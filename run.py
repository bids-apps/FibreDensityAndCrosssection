#!/usr/bin/env python3
import argparse
import os
import subprocess
from glob import glob
import sys
from lib.errorMessage  import errorMessage
from lib.printMessage  import printMessage
import lib.app, lib.cmdlineParser
from lib.runCommand    import runCommand
from lib.isWindows     import isWindows
from lib.getHeaderInfo import getHeaderInfo
from lib.getHeaderProperty import getHeaderProperty
from lib.delFile       import delFolder

__version__ = 'BIDS-App \'FibreDensityAndCrosssection\' version {}'.format(open('/version').read()) if os.path.exists('/version') else 'BIDS-App \'FibreDensityAndCrosssection\' standalone'


lib.app.author = 'David Raffelt (david.raffelt@florey.edu.au)'

lib.cmdlineParser.initialise('Perform group analysis of diffusion MRI data with a Fixel-Based Analysis (FBA) of Fibre Density, Fibre Cross-section and a combined measure (Fibre Density & Cross-section). The analysis pipeline relies primarily on the MRtrix3 software package (www.mrtrix.org).')

lib.app.parser.add_argument('bids_dir', help='The directory with the input dataset '
                                        'formatted according to the BIDS standard.')
lib.app.parser.add_argument('output_dir', help='The directory where the output files '
                                         'should be stored. If you are running group level analysis '
                                         'this folder should be prepopulated with the results of the '
                                         'participant level analysis.')

analysis_level_choices = ['participant1', 'group1', 'participant2', 'group2', 'participant3', 'group3', 'participant4', 'group4']

lib.app.parser.add_argument('analysis_level', help='Level of the analysis that will be performed. '
                                                   'Valid choices are: [' + ', '.join(analysis_level_choices) + ']. \nMultiple participant '
                                                   'level analyses can be run independently(in parallel) using the same output_dir.',
                                              choices = analysis_level_choices)

lib.app.parser.add_argument('--participant_label', help='The label(s) of the participant(s) that should be analyzed. The label '
                                                 'corresponds to sub-<participant_label> from the BIDS spec '
                                                 '(so it does not include "sub-"). If this parameter is not '
                                                 'provided all subjects should be analyzed. Multiple '
                                                 'participants can be specified with a space separated list.',
                                                  nargs='+')

lib.app.parser.add_argument('--n_cpus', type=int, default='1', help='The number of CPU cores available on the compute node. '
                                                             'Set to 0 to use the maximum number of cores available')

options = lib.app.parser.add_argument_group('Options for this Fibre Density and Cross-section BIDS-App')


options.add_argument('-vox_size', type=float, default='1.25', help='define the voxel size (in mm) to be used during the upsampling step (participant1 analysis level only)')

options.add_argument('-group_subset', help='Define a subset of participants to be used when generating the group-average FOD template and response functions. The subset is to be supplied as a comma separate list. Note the subset should be representable of your entire population and not biased towards one particular group. For example in a patient-control comparison, choose equal numbers of patients and controls. Used in group1 and group2 analysis levels.', nargs=1)

options.add_argument('-num_tracks', type=int, default='20000000', help='define the number of streamlines to be computed '
                                                                        'when performing tractography on the FOD template. '
                                                                        '(group3 analysis level only)')
options.add_argument('-num_tracks_sift', type=int, default='2000000', help='define the number of streamlines to '
                                                                           'remain after performing SIFT on the tractogram'
                                                                           '(group3 analysis level only)')



lib.app.initialise()

if isWindows():
  errorMessage('Script cannot be run on Windows due to FSL dependency')

nthreads = ''
if (lib.app.args.n_cpus):
  lib.app.mrtrixNThreads = '-nthreads ' + str(lib.app.args.n_cpus) + ' '

subjects_to_analyze = []
# only for a subset of subjects
if lib.app.args.participant_label:
  subjects_to_analyze = lib.app.args.participant_label
# for all subjects
else:
  subject_dirs = glob(os.path.join(lib.app.args.bids_dir, 'sub-*'))
  subjects_to_analyze = [subject_dir.split("-")[-1] for subject_dir in subject_dirs]

# create output subjects directory
all_subjects_dir = os.path.join(lib.app.args.output_dir, 'subjects');
if not os.path.exists(all_subjects_dir):
  os.mkdir(all_subjects_dir)

# create the output template directory
template_dir = os.path.join(lib.app.args.output_dir, 'template')
if not os.path.exists(template_dir):
  os.mkdir(template_dir)

# create a temporary directory for intermediate files
lib.app.makeTempDir()

# read in group subset if supplied
subset = []
if lib.app.args.group_subset:
  subset = lib.app.args.group_subset[0].split(',')


# running participant level 1 (basic preprocessing)
if lib.app.args.analysis_level == 'participant1':

  subprocess.check_call('bids-validator ' + lib.app.args.bids_dir, shell=True)

  for subject_label in subjects_to_analyze:
    label = 'sub-' + subject_label
    printMessage('running basic pre-processing for ' + label)

    # Read DWI(s) in BIDS folder
    all_dwi_images = glob(os.path.join(lib.app.args.bids_dir, label, '*dwi', '*_dwi.nii*'))

    # TODO handle multiple DWIs (e.g. time points) in subject directory
    if (len(all_dwi_images) > 1):
      errorMessage('Multiple DWIs found in subject folder. Multiple sessions not currently supported.')

    # Create output subject directory
    subject_dir = os.path.join(all_subjects_dir, subject_label)
    if not os.path.exists(subject_dir):
      os.mkdir(subject_dir)

    # Check existence output files from this analysis level
    dwi_preproc_file = os.path.join(subject_dir, 'dwi_preproc.mif')
    lib.app.checkOutputFile(dwi_preproc_file)
    wm_response_file = os.path.join(subject_dir, 'wm_response.txt')
    lib.app.checkOutputFile(wm_response_file)
    gm_response_file = os.path.join(subject_dir, 'gm_response.txt')
    lib.app.checkOutputFile(gm_response_file)
    csf_response_file = os.path.join(subject_dir, 'csf_response.txt')
    lib.app.checkOutputFile(csf_response_file)

    # DW gradient files
    grad_prefix = os.path.join(lib.app.args.bids_dir, label, 'dwi', label + '_dwi')
    if not (os.path.isfile(grad_prefix + '.bval') and os.path.isfile(grad_prefix + '.bvec')):
      grad_prefix = os.path.join(lib.app.args.bids_dir, 'dwi')
      if not (os.path.isfile(grad_prefix + '.bval') and os.path.isfile(grad_prefix + '.bvec')):
        errorMessage('Unable to locate valid diffusion gradient table');
    grad_import_option = ' -fslgrad ' + grad_prefix + '.bvec ' + grad_prefix + '.bval'

    # Json file
    json_path = os.path.join(lib.app.args.bids_dir, label, 'dwi', label + '_dwi.json')
    if os.path.isfile(json_path):
      json_import_option = ' -json_import ' + json_path
    else:
      json_import_option = ''


    # Stuff DWI gradients in *.mif file
    dwi_mrtrix_file = os.path.join(lib.app.tempDir, subject_label + 'dwi.mif')
    runCommand('mrconvert ' + all_dwi_images[0] + grad_import_option + json_import_option + ' ' + dwi_mrtrix_file)

    # Denoise
    dwi_denoised_file = os.path.join(lib.app.tempDir, subject_label + 'dwi_denoised.mif')
    runCommand('dwidenoise ' + dwi_mrtrix_file + ' ' + dwi_denoised_file)

    # Topup and eddy TODO add reverse phase encode capability
    runCommand('dwipreproc -rpe_none -pe_dir AP ' + lib.app.mrtrixNThreads + dwi_denoised_file + ' ' + dwi_preproc_file + lib.app.mrtrixForce)

    # Estimate WM, GM, CSF response functions
    runCommand('dwi2response ' + lib.app.mrtrixNThreads + ' dhollander ' + dwi_preproc_file + ' ' + wm_response_file + ' '
                                                                                                  + gm_response_file + ' '
                                                                                                  + csf_response_file + lib.app.mrtrixForce)

# running group level 1 (average response functions) TODO check for user supplied subset to ensure response is not biased
elif lib.app.args.analysis_level == "group1":
  printMessage('averaging response functions')

  # Check output files exist
  wm_response_file = os.path.join(lib.app.args.output_dir, 'average_wm_response.txt')
  lib.app.checkOutputFile(wm_response_file)
  gm_response_file = os.path.join(lib.app.args.output_dir, 'average_gm_response.txt')
  lib.app.checkOutputFile(gm_response_file)
  csf_response_file = os.path.join(lib.app.args.output_dir, 'average_csf_response.txt')
  lib.app.checkOutputFile(csf_response_file)

  # process subset
  input_wm_files = []
  input_gm_files = []
  input_csf_files = []
  if (len(subset) > 0):
    printMessage('Using a group subset to compute average response functions' + str(subset))
    subject_labels = [os.path.basename(x) for x in glob(os.path.join(all_subjects_dir, '*'))]
    for subj in subset:
      if subj not in subject_labels:
        errorMessage('subject label (' + os.path.basename(subj) + ') supplied as part of -group_subset option does exist in subjects directory')
      input_wm_files.append(os.path.join(all_subjects_dir, subj, 'wm_response.txt'))
      input_gm_files.append(os.path.join(all_subjects_dir, subj, 'gm_response.txt'))
      input_csf_files.append(os.path.join(all_subjects_dir, subj, 'csf_response.txt'))
  # use all subjects
  else:
    input_wm_files = glob(os.path.join(all_subjects_dir, '*', 'wm_response.txt'))
    input_gm_files = glob(os.path.join(all_subjects_dir, '*', 'gm_response.txt'))
    input_csf_files = glob(os.path.join(all_subjects_dir, '*', 'csf_response.txt'))

  runCommand('average_response ' + ' '.join(input_wm_files) + ' ' + wm_response_file + lib.app.mrtrixForce)
  runCommand('average_response ' + ' '.join(input_gm_files) + ' ' + gm_response_file + lib.app.mrtrixForce)
  runCommand('average_response ' + ' '.join(input_csf_files) + ' ' + csf_response_file + lib.app.mrtrixForce)


# running participant level 2 (upsample, compute brain masks and FODs, perform intensity normalisation and bias field correction)
elif lib.app.args.analysis_level == "participant2":
  for subject_label in subjects_to_analyze:

    subject_dir = os.path.join(all_subjects_dir, subject_label)
    output_mask = os.path.join(subject_dir, 'mask.mif')
    lib.app.checkOutputFile(output_mask)
    output_fod = os.path.join(subject_dir, 'fod.mif')
    lib.app.checkOutputFile(output_fod)
    output_gm = os.path.join(subject_dir, 'gm.mif')
    lib.app.checkOutputFile(output_gm)
    output_csf = os.path.join(subject_dir, 'csf.mif')
    lib.app.checkOutputFile(output_csf)

    min_voxel_size = 1.25;
    if lib.app.args.vox_size:
      min_voxel_size = float(lib.app.args.vox_size)

    voxel_sizes = getHeaderInfo(os.path.join(subject_dir, 'dwi_preproc.mif'), 'vox').split()
    mean_voxel_size = 0.0
    for i in range(0,3):
      mean_voxel_size = mean_voxel_size + float(voxel_sizes[i]) / 3.0

    input_to_csd = os.path.join(subject_dir, 'dwi_preproc.mif')

    # Up sample
    if mean_voxel_size > min_voxel_size:
      runCommand('mrresize -vox ' + str(min_voxel_size) + ' ' + input_to_csd + ' ' + os.path.join(lib.app.tempDir, subject_label + 'dwi_upsampled.mif'))
      input_to_csd = os.path.join(lib.app.tempDir, subject_label + 'dwi_upsampled.mif')

    # Compute brain mask
    runCommand('dwi2mask ' + input_to_csd + ' ' + output_mask + lib.app.mrtrixForce)

    # Perform CSD
    runCommand('dwi2fod msmt_csd ' + input_to_csd + ' -mask ' + output_mask + ' ' +
                os.path.join(lib.app.args.output_dir, 'average_wm_response.txt') + ' ' +  os.path.join(lib.app.tempDir, subject_label + 'fod.mif') + ' ' +
                os.path.join(lib.app.args.output_dir, 'average_gm_response.txt') + ' ' + os.path.join(lib.app.tempDir, subject_label + 'gm.mif') + ' ' +
                os.path.join(lib.app.args.output_dir, 'average_csf_response.txt') + ' ' + os.path.join(lib.app.tempDir, subject_label + 'csf.mif'))

    runCommand('mtbin -independent ' + os.path.join(lib.app.tempDir, subject_label + 'fod.mif') + ' ' + output_fod + ' ' +
                                       os.path.join(lib.app.tempDir, subject_label + 'gm.mif')  + ' ' + output_gm  + ' ' +
                                       os.path.join(lib.app.tempDir, subject_label + 'csf.mif') + ' ' + output_csf + lib.app.mrtrixForce)


# running group level 2 (generate FOD template)
elif lib.app.args.analysis_level == 'group2':

  # TODO if user supplies a subset, then only output the template

  # Check if outputs exist
  fod_template = os.path.join(template_dir, 'fod_template.mif')
  lib.app.checkOutputFile(fod_template)

  lib.app.gotoTempDir()
  os.mkdir('fod_input')
  os.mkdir('mask_input')

  # Check if all members of subset exist
  if (len(subset) > 0):
    printMessage('Using a group subset to compute population template' + str(subset))
    subject_labels = [os.path.basename(x) for x in glob(os.path.join(all_subjects_dir, '*'))]
    for subj in subset:
      if subj not in subject_labels:
        errorMessage('subject label (' + os.path.basename(subj) + ') supplied as part of -group_subset option does exist in subjects directory')
      lib.app.checkOutputFile(os.path.join(all_subjects_dir, subj, 'subject2template_warp.mif'))
      lib.app.checkOutputFile(os.path.join(all_subjects_dir, subj, 'template2subject_warp.mif'))


  # make symlinks to all population_template inputs in single directory
  for subj in glob(os.path.join(all_subjects_dir, '*')):
    if (len(subset) > 0):
      if os.path.basename(subj) in subset:
        os.symlink(os.path.join(subj, 'fod.mif'), os.path.join('fod_input', os.path.basename(subj) + '.mif'))
        os.symlink(os.path.join(subj, 'mask.mif'), os.path.join('mask_input', os.path.basename(subj) + '.mif'))
    else:
      os.symlink(os.path.join(subj, 'fod.mif'), os.path.join('fod_input', os.path.basename(subj) + '.mif'))
      os.symlink(os.path.join(subj, 'mask.mif'), os.path.join('mask_input', os.path.basename(subj) + '.mif'))

  # Compute FOD template
  if (len(subset) > 0):
    runCommand('population_template fod_input -mask mask_input ' + os.path.join(lib.app.tempDir, 'tmp.mif'))
    # Set a field in the header of the template to mark it as being generated as a subset or not. This is used in the next step
    runCommand('mrconvert ' + os.path.join(lib.app.tempDir, 'tmp.mif') + ' -header_set made_from_subset true ' + fod_template + lib.app.mrtrixForce)
  else:
    runCommand('population_template fod_input -mask mask_input ' + os.path.join(lib.app.tempDir, 'tmp.mif') + ' -warp_dir ' +  os.path.join(lib.app.tempDir, 'warps'))
    runCommand('mrconvert ' + os.path.join(lib.app.tempDir, 'tmp.mif') + ' -header_set made_from_subset false ' + fod_template + lib.app.mrtrixForce)
    # Save all warps since we don't need to generate them in the next step if all subjects were used to make the template
    for subj in [os.path.basename(x) for x in glob(os.path.join(all_subjects_dir, '*'))]:
      runCommand('warpconvert -type warpfull2deformation -template ' + fod_template + ' '
                              + os.path.join(lib.app.tempDir, 'warps', subj + '.mif') + ' '
                              + os.path.join(all_subjects_dir, subj, 'subject2template_warp.mif') + lib.app.mrtrixForce)
      runCommand('warpconvert -type warpfull2deformation -from 2 -template ' + os.path.join(all_subjects_dir, subj, 'fod.mif') + ' '
                              + os.path.join(lib.app.tempDir, 'warps', subj + '.mif') + ' '
                              + os.path.join(all_subjects_dir, subj, 'template2subject_warp.mif') + lib.app.mrtrixForce)


# running participant level 3 (register FODs, and warp masks)
elif lib.app.args.analysis_level == "participant3":

  for subject_label in subjects_to_analyze:
    subject_dir = os.path.join(all_subjects_dir, subject_label)

    # Check existence of output
    mask_template = os.path.join(subject_dir, 'mask_in_template_space.mif')
    lib.app.checkOutputFile(mask_template)

    # if the template was generated with a subset of the whole study, then we still need to register all subjects to that template
    if getHeaderProperty (os.path.join(template_dir, 'fod_template.mif'), 'made_from_subset') == 'true':
      subject2template = os.path.join(subject_dir, 'subject2template_warp.mif')
      lib.app.checkOutputFile(subject2template)
      template2subject = os.path.join(subject_dir, 'template2subject_warp.mif')
      lib.app.checkOutputFile(template2subject)
      # TODO If only a subset of images were used to make the template, then register all images to the template
      runCommand('mrregister ' + os.path.join(subject_dir, 'fod.mif') + ' -mask1 ' + os.path.join(subject_dir, 'mask.mif')
                               + ' ' +  os.path.join(template_dir, 'fod_template.mif') + ' -nl_warp '
                               + subject2template + ' ' + template2subject + lib.app.mrtrixForce)


    # Transform masks into template space. This is used in the group3 analysis level for trimming the
    # final voxek mask to exclude voxels that do not contain data from all subjects
    runCommand('mrtransform ' + os.path.join(subject_dir, 'mask.mif') + ' -warp ' + subject2template + ' -interp nearest '
                              + mask_template + lib.app.mrtrixForce)


# running group level 3 compute voxel and fixel masks, tractography, sift)
elif lib.app.args.analysis_level == "group3":

  voxel_mask = os.path.join(template_dir, 'voxel_mask.mif')
  lib.app.checkOutputFile(voxel_mask)
  fixel_mask = os.path.join(template_dir, 'fixel_mask')
  lib.app.checkOutputFile(fixel_mask)
  tracks = os.path.join(template_dir, 'tracks.tck')
  lib.app.checkOutputFile(tracks)
  tracks_sift = os.path.join(template_dir, 'tracks_sift.tck')
  lib.app.checkOutputFile(tracks_sift)
  fod_template = os.path.join(template_dir, 'fod_template.mif')

  # Check tractography and SIFT input
  num_tracks_sift = 2000000;
  if lib.app.args.num_tracks:
    num_tracks_sift = int(lib.app.args.num_tracks_sift)
  num_tracks = 20000000;
  if lib.app.args.num_tracks:
    num_tracks = int(lib.app.args.num_tracks)
  if num_tracks_sift >= num_tracks:
    errorMessage('the tracks remaining after SIFT must be less than the number of tracks generated during tractography')

  # Compute voxel mask and intersect with all brain masks to ensure mask voxels are present in all subjects
  all_masks_in_template_space = glob(os.path.join(all_subjects_dir, '*', 'mask_in_template_space.mif'))
  runCommand('mrconvert -coord 3 0 ' + fod_template + ' - | mrthreshold - - | mrmath - ' + ' '.join(all_masks_in_template_space)
             + ' min ' + voxel_mask + lib.app.mrtrixForce)

  # Compute fixel mask
  delFolder(fixel_mask)
  runCommand('fod2fixel -mask ' + voxel_mask + ' -fmls_peak_value 0.2 ' + fod_template + ' ' + fixel_mask + lib.app.mrtrixForce)

  # Perform tractography on the FOD template
  runCommand('tckgen -angle 22.5 -maxlen 250 -minlen 10 -power 1.0 ' + fod_template + ' -seed_image '
             + voxel_mask + ' -mask ' + voxel_mask + ' -number ' + str(num_tracks) + ' ' + tracks + lib.app.mrtrixForce)

  # SIFT the streamlines
  runCommand('tcksift ' + tracks + ' ' + fod_template + ' -term_number ' + str(num_tracks_sift) + ' ' + tracks_sift + lib.app.mrtrixForce)


# Warp FODs, Compute FD, reorient fixels, fixelcorrespondence , compute FC, compute FDC per subject
elif lib.app.args.analysis_level == "participant4":

  for subject_label in subjects_to_analyze:
    subject_dir = os.path.join(all_subjects_dir, subject_label)

    # Fixel output all aligned in template space
    template_fd = os.path.join(template_dir, 'fd')
    lib.app.checkOutputFile(template_fd)
    template_fc = os.path.join(template_dir, 'fc')
    lib.app.checkOutputFile(template_fc)
    template_log_fc = os.path.join(template_dir, 'log_fc', subject_label + '.mif')
    lib.app.checkOutputFile(template_log_fc)
    template_fdc = os.path.join(template_dir, 'fdc', subject_label + '.mif')
    lib.app.checkOutputFile(template_fdc)


    # Transform FOD images (without reorientation)
    runCommand('mrtransform -noreorientation ' + os.path.join(subject_dir, 'fod.mif')
                                               + ' -warp ' + os.path.join(subject_dir, 'subject2template_warp.mif') + ' '
                                               + os.path.join(lib.app.tempDir, subject_label + 'fod_warped.mif'))


    # Segment each FOD into fixels and compute AFD integral per fixel
    delFolder(os.path.join(lib.app.tempDir, subject_label + 'fixel_fd'))
    runCommand('fod2fixel -afd fd.mif ' + os.path.join(lib.app.tempDir, subject_label + 'fod_warped.mif')
                                        + ' -mask ' + os.path.join(template_dir, 'voxel_mask.mif') + ' '
                                        + os.path.join(lib.app.tempDir, subject_label + 'fixel_fd'))

    # Reorient each fixel's direction (inplace) to account for the spatial tranformation
    runCommand('fixelreorient -force ' + os.path.join(lib.app.tempDir, subject_label + 'fixel_fd') + ' '
                                       + os.path.join(subject_dir, 'subject2template_warp.mif') + ' '
                                       + os.path.join(lib.app.tempDir, subject_label + 'fixel_fd'))


    # For each fixel in template space, find the corresponding fixel in the subject and assign its AFD value.
    # Put all subjects AFD in the sample folder under the template directory ready for statistical analysis
    runCommand('fixelcorrespondence ' + os.path.join(lib.app.tempDir, subject_label + 'fixel_fd', 'fd.mif') + ' '
                                      + os.path.join(template_dir, 'fixel_mask') + ' '
                                      + template_fd + ' ' + subject_label + '.mif' + lib.app.mrtrixForce)

    # Compute FC
    runCommand('warp2metric ' + os.path.join(subject_dir, 'subject2template_warp.mif')
                              + ' -fc ' + os.path.join(template_dir, 'fixel_mask') + ' '
                              + template_fc + ' ' + subject_label + '.mif' + lib.app.mrtrixForce)

    # Compute log FC
    if not os.path.exists(os.path.join(template_dir, 'log_fc')):
      os.mkdir(os.path.join(template_dir, 'log_fc'))
    runCommand('mrcalc ' + os.path.join(template_fc, subject_label + '.mif') + ' -log ' + template_log_fc + lib.app.mrtrixForce)
    if not os.path.exists(os.path.join(template_dir, 'fdc')):
      os.mkdir(os.path.join(template_dir, 'fdc'))
    runCommand('mrcalc ' + os.path.join(template_fd, subject_label + '.mif') + ' ' + os.path.join(template_fc, subject_label + '.mif')
                         + ' -mult ' + template_fdc + lib.app.mrtrixForce)

  # Copy index and directions file into log FC and FDC fixel directories
  runCommand('cp ' + os.path.join(template_fc, 'index.mif') + ' '
                   + os.path.join(template_fc, 'directions.mif') + ' '
                   + os.path.join(template_dir, 'log_fc'))
  runCommand('cp ' + os.path.join(template_fc, 'index.mif') + ' '
                   + os.path.join(template_fc, 'directions.mif') + ' '
                   + os.path.join(template_dir, 'fdc'))

# Perform fixel-based statistical inference in FD, FC and FDC
elif lib.app.args.analysis_level == "group4":
  print('asdf')


lib.app.complete()
