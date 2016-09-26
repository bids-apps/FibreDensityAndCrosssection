#!/usr/bin/env python3
import argparse
import os
import subprocess
from glob import glob
import sys
from lib.errorMessage  import errorMessage
import lib.app, lib.cmdlineParser
from lib.runCommand    import runCommand

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
options.add_argument('--n_cpus', type=int, default='1', help='The number of CPU cores available on the compute node. '
                                                             'Set to 0 to use the maximum number of cores available')

lib.app.initialise()

subprocess.check_call('bids-validator ' + lib.app.args.bids_dir, shell=True)

nthreads = ""
if (lib.app.args.n_cpus):
  nthreads = "-nthreads " + str(lib.app.args.n_cpus)

subjects_to_analyze = []
# only for a subset of subjects
if lib.app.args.participant_label:
  subjects_to_analyze = lib.app.args.participant_label
# for all subjects
else:
  subject_dirs = glob(os.path.join(lib.app.args.bids_dir, "sub-*"))
  subjects_to_analyze = [subject_dir.split("-")[-1] for subject_dir in subject_dirs]

# running participant level 1 (basic preprocessing)
if lib.app.args.analysis_level == "participant1":

  # dwipreproc
  for subject_label in subjects_to_analyze:
    label = 'sub-' + subject_label
    all_dwi_images = glob(os.path.join(lib.app.args.bids_dir, label, "*dwi", "*_dwi.nii*"))

    # TODO handle multiple DWIs in subject folder
    if (len(all_dwi_images) > 1):
      errorMessage("Multiple DWIs found in subject folder. Multiple sessions not currently supported.")

    if not os.path.exists(os.path.join(lib.app.args.output_dir, "dwi")):
      os.mkdir(os.path.join(lib.app.args.output_dir, "dwi"))

    grad_prefix = os.path.join(lib.app.args.bids_dir, label, 'dwi', label + '_dwi')
    if not (os.path.isfile(grad_prefix + '.bval') and os.path.isfile(grad_prefix + '.bvec')):
      grad_prefix = os.path.join(lib.app.args.bids_dir, 'dwi')
      if not (os.path.isfile(grad_prefix + '.bval') and os.path.isfile(grad_prefix + '.bvec')):
        errorMessage('Unable to locate valid diffusion gradient table');
    grad_import_option = ' -fslgrad ' + grad_prefix + '.bvec ' + grad_prefix + '.bval'
    json_path = os.path.join(lib.app.args.bids_dir, label, 'dwi', label + '_dwi.json')
    if os.path.isfile(json_path):
      json_import_option = ' -json_import ' + json_path
    else:
      json_import_option = ''

    # Stuff gradients in mif file
    dwi_mrtrix_file = os.path.join(lib.app.args.output_dir, 'dwi', subject_label + '_input.mif')
    runCommand('mrconvert ' + all_dwi_images[0] + grad_import_option + json_import_option + ' ' + dwi_mrtrix_file)

    # Denoise
    dwi_denoised_file = os.path.join(lib.app.args.output_dir, 'dwi', subject_label + '_denoised.mif')
    runCommand('dwidenoise ' + dwi_mrtrix_file + ' ' + dwi_denoised_file)

    # topup and eddy
    dwi_preproc_file = os.path.join(lib.app.args.output_dir, 'dwi', subject_label + '_preproc.mif')
    #runCommand("dwipreproc %s -rpe_none AP %s %s"%(nthreads, 'dwi_denoised.mif', dwi_preproc_file))

    # Compute brain mask
    if not os.path.exists(os.path.join(lib.app.args.output_dir, "mask")):
      os.mkdir(os.path.join(lib.app.args.output_dir, "mask"))
    mask_file = os.path.join(lib.app.args.output_dir, "mask", subject_label + '_mask.mif')
    #runCommand("dwi2mask -nthreads %s %s %s"%(nthreads, preproc_file, mask_file))

    # Perform bias field correction
    bias_file = os.path.join(lib.app.args.output_dir, 'dwi', subject_label + '_preproc_bias.mif')
    #runCommand("dwibiascorrect -ants -mask %s %s %s"%(mask_file, preproc_file, bias_file))

# running group level 1 (intensity normalisation)
elif lib.app.args.analysis_level == "group1":
  mask_dir = os.path.join(lib.app.args.output_dir, "mask")
  dwi_dir = os.path.join(lib.app.args.output_dir, "dwi")
  norm_dwi_dir = os.path.join(lib.app.args.output_dir, "norm_dwi")
  fa_template = os.path.join(lib.app.args.output_dir, "fa_template.mif")
  wm_mask = os.path.join(lib.app.args.output_dir, "wm_mask.mif")
  cmd = "dwiintensitynorm %s %s %s %s %s %s"%(nthreads, dwi_dir, mask_dir, norm_dwi_dir, fa_template, wm_mask)
  print (cmd)
  subprocess.check_call(cmd, shell=True)
  response_dir = os.path.join(lib.app.args.output_dir, "response")
  if not os.path.exists(response_dir):
    os.mkdir(response_dir)

# running participant level 2 (dwi2response)
elif lib.app.args.analysis_level == "participant2":
  for subject_label in subjects_to_analyze:
    input_dwi = os.path.join(lib.app.args.output_dir, "norm_dwi", "sub-" + subject_label + "_dwi.mif")
    output_response = os.path.join(lib.app.args.output_dir, "response", subject_label + "-response.txt")
    cmd = "dwi2response %s tournier %s %s"%(nthreads, input_dwi, output_response)
    print (cmd)
    subprocess.check_call(cmd, shell=True)

# running group level 2 (average response)
elif lib.app.args.analysis_level == "group2":
  cmd = "average_response %s/* %s"%(os.path.join(lib.app.args.output_dir, "response"), os.path.join(lib.app.args.output_dir, "average_response.txt"))
  subprocess.check_call(cmd, shell=True)

# running participant level 3 (upsample, brain mask, CSD)
elif lib.app.args.analysis_level == "participant3":
  for subject_label in subjects_to_analyze:
    input_dwi = os.path.join(lib.app.args.output_dir, "norm_dwi", "sub-" + subject_label + "_dwi.mif")
    if not os.path.exists(os.path.join(lib.app.args.output_dir, "upsampled_dwi")):
      os.mkdir(os.path.join(lib.app.args.output_dir, "upsampled_dwi"))
    upsampled_dwi = os.path.join(lib.app.args.output_dir, "upsampled_dwi", "sub-" + subject_label + "_dwi.mif")
    #TODO, detect if resolution > 2mm
    cmd = "mrresize %s -scale 2.0 %s %s"%(nthreads, input_dwi, upsampled_dwi)
    subprocess.check_call(cmd, shell=True)
    if not os.path.exists(os.path.join(lib.app.args.output_dir, "upsampled_mask")):
      os.mkdir(os.path.join(lib.app.args.output_dir, "upsampled_mask"))
    input_mask = os.path.join(lib.app.args.output_dir, "mask", "sub-" + subject_label + "_mask.mif")
    upsampled_mask = os.path.join(lib.app.args.output_dir, "upsampled_mask", "sub-" + subject_label + "_mask.mif")
    cmd = "mrresize %s -scale 2.0 -interp nearest %s %s"%(nthreads, input_mask, upsampled_mask)
    subprocess.check_call(cmd, shell=True)
    if not os.path.exists(os.path.join(lib.app.args.output_dir, "fod")):
      os.mkdir(os.path.join(lib.app.args.output_dir, "fod"))
    fod = os.path.join(lib.app.args.output_dir, "fod", "sub-" + subject_label + "_fod.mif")
    group_response = os.path.join(lib.app.args.output_dir, "response", subject_label + "-response.txt")
    cmd = "dwiextract %s %s - | dwi2fod %s msmt_csd - %s %s -mask %s"%(nthreads, upsampled_dwi, nthreads, group_response, fod, upsampled_mask)
    subprocess.check_call(cmd, shell=True)

# running group level 3 (generate study specific FOD template, compute fixel template mask, compute voxel mask, tractography, SIFT)
elif lib.app.args.analysis_level == "group3":
  # TODO if more than 20 subjects randomly select subset to generate template
  fod = os.path.join(lib.app.args.output_dir, "fod")
  upsampled_mask = os.path.join(lib.app.args.output_dir, "upsampled_mask")
  fod_template = os.path.join(lib.app.args.output_dir, "fod_template.mif")
  # TODO remove level adjustment after test
  cmd = "population_template %s -nl_scale 0.5,0.75,1.0 -nl_lmax 2,2,2 -nl_niter 5,5,5 %s -mask_dir %s %s"%(nthreads, fod, upsampled_mask, fod_template)
  subprocess.check_call(cmd, shell=True)

# running participant level 3 (register FODs, warp FODs, compute FD, reorient fixels, fixelcorrespondence, compute FC, compute FDC)
elif lib.app.args.analysis_level == "participant4":
  if not os.path.exists(os.path.join(lib.app.args.output_dir, "warp")):
    os.mkdir(os.path.join(lib.app.args.output_dir, "warp"))
  if not os.path.exists(os.path.join(lib.app.args.output_dir, "warped_fod")):
      os.mkdir(os.path.join(lib.app.args.output_dir, "warped_fod"))
  fod_template = os.path.join(lib.app.args.output_dir, "fod_template.mif")
  for subject_label in subjects_to_analyze:
    fod = os.path.join(lib.app.args.output_dir, "fod", "sub-" + subject_label + "_fod.mif")
    upsampled_mask = os.path.join(lib.app.args.output_dir, "upsampled_mask", "sub-" + subject_label + "_mask.mif")
    subject2template_warp = os.path.join(lib.app.args.output_dir, "warp", "sub-" + subject_label + "subject2template_warp.mif")
    template2subject_warp = os.path.join(lib.app.args.output_dir, "warp", "sub-" + subject_label + "template2subject_warp.mif")
    cmd = "mrregister %s %s -mask1 %s %s -nl_warp %s %s"%(nthreads, fod, upsampled_mask, fod_template, subject2template_warp, template2subject_warp)
    subprocess.check_call(cmd, shell=True)
    warped_fod = os.path.join(lib.app.args.output_dir, "warped_fod", "sub-" + subject_label + "_warped_fod.mif")
    cmd = "mrtransform %s -noreorientation %s -warp %s %s"%(nthreads, fod, subject2template_warp, warped_fod)
    subprocess.check_call(cmd, shell=True)


# running group level 3 (fixelcfestats)
elif lib.app.args.analysis_level == "group4":
  print ("group4")
