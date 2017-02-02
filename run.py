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


options = lib.app.parser.add_argument_group('Options for the population_template script')


options.add_argument('--participant_label', help='The label(s) of the participant(s) that should be analyzed. The label '
                                                 'corresponds to sub-<participant_label> from the BIDS spec '
                                                 '(so it does not include "sub-"). If this parameter is not '
                                                 'provided all subjects should be analyzed. Multiple '
                                                 'participants can be specified with a space separated list.',
                                                  nargs=1)
options.add_argument('--n_cpus', type=int, default='+', help='The number of CPU cores available on the compute node. '
                                                             'Set to 0 to use the maximum number of cores available')

options.add_argument('--vox_size', type=float, default='1.25', help='define the voxel size (in mm) to be used during the upsampling step')




lib.app.initialise()

if isWindows():
  errorMessage('Script cannot be run on Windows due to FSL dependency')

subprocess.check_call('bids-validator ' + lib.app.args.bids_dir, shell=True)

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

# running participant level 1 (basic preprocessing)
if lib.app.args.analysis_level == 'participant1':

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

    lib.app.makeTempDir()
    # Stuff DWI gradients in *.mif file
    dwi_mrtrix_file = os.path.join(lib.app.tempDir, 'dwi.mif')
    runCommand('mrconvert ' + all_dwi_images[0] + grad_import_option + json_import_option + ' ' + dwi_mrtrix_file)

    # Denoise
    dwi_denoised_file = os.path.join(lib.app.tempDir, 'dwi_denoised.mif')
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
  wm_response_file = os.path.join(all_subjects_dir, 'average_wm_response.txt')
  lib.app.checkOutputFile(wm_response_file)
  gm_response_file = os.path.join(all_subjects_dir, 'average_gm_response.txt')
  lib.app.checkOutputFile(gm_response_file)
  csf_response_file = os.path.join(all_subjects_dir, 'average_csf_response.txt')
  lib.app.checkOutputFile(csf_response_file)
  input_wm_files = glob(os.path.join(all_subjects_dir, '*', 'wm_response.txt'))
  runCommand('average_response ' + ' '.join(input_wm_files) + ' ' + wm_response_file + lib.app.mrtrixForce)
  input_gm_files = glob(os.path.join(all_subjects_dir, '*', 'gm_response.txt'))
  runCommand('average_response ' + ' '.join(input_gm_files) + ' ' + gm_response_file + lib.app.mrtrixForce)
  input_csf_files = glob(os.path.join(all_subjects_dir, '*', 'csf_response.txt'))
  runCommand('average_response ' + ' '.join(input_csf_files) + ' ' + csf_response_file + lib.app.mrtrixForce)


# running participant level 2 (upsample, compute brain masks and FODs)
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
      lib.app.makeTempDir()
      runCommand('mrresize -vox ' + str(min_voxel_size) + ' ' + input_to_csd + ' ' + os.path.join(lib.app.tempDir, 'dwi_upsampled.mif'))
      input_to_csd = os.path.join(lib.app.tempDir, 'dwi_upsampled.mif')

    # Compute brain mask
    runCommand('dwi2mask ' + input_to_csd + ' ' + output_mask + lib.app.mrtrixForce)

    # Perform CSD
    runCommand('dwi2fod msmt_csd ' + input_to_csd + ' -mask ' + output_mask + ' ' +
                os.path.join(all_subjects_dir, 'average_wm_response.txt') + ' ' +  os.path.join(lib.app.tempDir, 'fod.mif') + ' ' +
                os.path.join(all_subjects_dir, 'average_gm_response.txt') + ' ' + os.path.join(lib.app.tempDir, 'gm.mif') + ' ' +
                os.path.join(all_subjects_dir, 'average_csf_response.txt') + ' ' + os.path.join(lib.app.tempDir, 'csf.mif'))

    runCommand('mtbin -independent ' + os.path.join(lib.app.tempDir, 'fod.mif') + ' ' + output_fod + ' ' +
                                       os.path.join(lib.app.tempDir, 'gm.mif')  + ' ' + output_gm  + ' ' +
                                       os.path.join(lib.app.tempDir, 'csf.mif') + ' ' + output_csf + lib.app.mrtrixForce)


# running group level 2 (generate FOD template, compute voxel and fixel masks, tractography and SIFT)
elif lib.app.args.analysis_level == 'group2':

  # TODO if user supplies a subset, then only output the template
  template_dir = os.path.join(all_subjects_dir, 'template')
  if not os.path.exists(template_dir):
    os.mkdir(template_dir)

  # Check if outputs exist
  fod_template = os.path.join(template_dir, 'fod_template.mif')
  lib.app.checkOutputFile(fod_template)
  voxel_mask = os.path.join(template_dir, 'voxel_mask.mif')
  lib.app.checkOutputFile(voxel_mask)
  fixel_mask = os.path.join(template_dir, 'fixel_mask')
  lib.app.checkOutputFile(fixel_mask)
  tracks = os.path.join(template_dir, 'tracks_20_million.tck')
  lib.app.checkOutputFile(tracks)
  tracks_sift = os.path.join(template_dir, 'tracks_2_million_sift.tck')
  lib.app.checkOutputFile(tracks_sift)

  lib.app.makeTempDir()
  lib.app.gotoTempDir()
  os.mkdir('fod_input')
  os.mkdir('mask_input')

  # make symlinks to all population_template inputs in single directory
  for subj in glob(os.path.join(all_subjects_dir, '*')):
    print(os.path.basename(subj))
    os.symlink(os.path.join(subj, 'fod.mif'), os.path.join('fod_input', os.path.basename(subj) + '.mif'))
    os.symlink(os.path.join(subj, 'mask.mif'), os.path.join('mask_input', os.path.basename(subj) + '.mif'))

  # Compute FOD template
  runCommand('population_template fod_input -mask mask_input ' + fod_template + lib.app.mrtrixForce)

  # Compute voxel mask
  runCommand('mrconvert -coord 3 0 ' + fod_template + ' - | mrthreshold - ' + voxel_mask + lib.app.mrtrixForce)

  # Compute fixel mask
  runCommand('fod2fixel -mask ' + voxel_mask + ' -fmls_peak_value 0.2 ' + fod_template + ' ' + fixel_mask + lib.app.mrtrixForce)

# running participant level 3 (register FODs, warp FODs, compute FD, reorient fixels, fixelcorrespondence, compute FC, compute FDC)
#elif lib.app.args.analysis_level == "participant3":
#  for subject_label in subjects_to_analyze:
#    subject_dir = os.path.join(all_subjects_dir, subject_label)


# running group level 3 (fixelcfestats
#elif lib.app.args.analysis_level == "group3":
#  # TODO if more than 20 subjects randomly select subset to generate template
#  fod = os.path.join(lib.app.args.output_dir, "fod")
#  upsampled_mask = os.path.join(lib.app.args.output_dir, "upsampled_mask")
#  fod_template = os.path.join(lib.app.args.output_dir, "fod_template.mif")
#  # TODO remove level adjustment after test
#  cmd = "population_template %s -nl_scale 0.5,0.75,1.0 -nl_lmax 2,2,2 -nl_niter 5,5,5 %s -mask_dir %s %s"%(nthreads, fod, upsampled_mask, fod_template)
#  subprocess.check_call(cmd, shell=True)

# running participant level 3 (register FODs, warp FODs, compute FD, reorient fixels, fixelcorrespondence, compute FC, compute FDC)
#elif lib.app.args.analysis_level == "participant4":
#  if not os.path.exists(os.path.join(lib.app.args.output_dir, "warp")):
#    os.mkdir(os.path.join(lib.app.args.output_dir, "warp"))
#  if not os.path.exists(os.path.join(lib.app.args.output_dir, "warped_fod")):
#      os.mkdir(os.path.join(lib.app.args.output_dir, "warped_fod"))
#  fod_template = os.path.join(lib.app.args.output_dir, "fod_template.mif")
#  for subject_label in subjects_to_analyze:
#    fod = os.path.join(lib.app.args.output_dir, "fod", "sub-" + subject_label + "_fod.mif")
#    upsampled_mask = os.path.join(lib.app.args.output_dir, "upsampled_mask", "sub-" + subject_label + "_mask.mif")
#    subject2template_warp = os.path.join(lib.app.args.output_dir, "warp", "sub-" + subject_label + "subject2template_warp.mif")
#    template2subject_warp = os.path.join(lib.app.args.output_dir, "warp", "sub-" + subject_label + "template2subject_warp.mif")
#    cmd = "mrregister %s %s -mask1 %s %s -nl_warp %s %s"%(nthreads, fod, upsampled_mask, fod_template, subject2template_warp, template2subject_warp)
#    subprocess.check_call(cmd, shell=True)
#    warped_fod = os.path.join(lib.app.args.output_dir, "warped_fod", "sub-" + subject_label + "_warped_fod.mif")
#    cmd = "mrtransform %s -noreorientation %s -warp %s %s"%(nthreads, fod, subject2template_warp, warped_fod)
#    subprocess.check_call(cmd, shell=True)


lib.app.complete()
