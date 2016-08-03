#!/usr/bin/env python3
import argparse
import os
import subprocess
from glob import glob
import sys

parser = argparse.ArgumentParser(description='Fixel-based analysis of fibre density (FD), fibre cross-section (FC), and fibre density and cross-section (FDC).')
parser.add_argument('bids_dir', help='The directory with the input dataset '
                                     'formatted according to the BIDS standard.')
parser.add_argument('output_dir', help='The directory where the output files '
                                       'should be stored. If you are running group level analysis '
                                       'this folder should be prepopulated with the results of the'
                                       'participant level analysis.')
parser.add_argument('analysis_level', help='Level of the analysis that will be performed. '
                                      'Multiple participant level analyses can be run independently '
                                      '(in parallel) using the same output_dir.',
                                      choices=['participant1', 'group1', 'participant2', 'group2', 'participant3', 'group3', 'participant4', 'group4'])
parser.add_argument('--participant_label', help='The label(s) of the participant(s) that should be analyzed. The label '
                                           'corresponds to sub-<participant_label> from the BIDS spec '
                                           '(so it does not include "sub-"). If this parameter is not '
                                           'provided all subjects should be analyzed. Multiple '
                                           'participants can be specified with a space separated list.',
                                            nargs="+")
parser.add_argument('--n_cpus', type=int, default='2', help='The number of CPUs available on the compute node')  #TODO pass to each cmd

args = parser.parse_args()

subjects_to_analyze = []
# only for a subset of subjects
if args.participant_label:
  subjects_to_analyze = args.participant_label
# for all subjects
else:
  subject_dirs = glob(os.path.join(args.bids_dir, "sub-*"))
  subjects_to_analyze = [subject_dir.split("-")[-1] for subject_dir in subject_dirs]

# running participant level 1 (basic preprocessing)
if args.analysis_level == "participant1":

  # TODO add detection of reverse phase encode b0s
  # dwipreproc
  for subject_label in subjects_to_analyze:
    for dwi_file in glob(os.path.join(args.bids_dir, "sub-%s"%subject_label,
                     "dwi", "*_dwi.nii*")) + glob(os.path.join(args.bids_dir,"sub-%s"%subject_label,"ses-*","dwi", "*dwi.nii*")):
      bvec = os.path.join(args.bids_dir, "sub-%s"%subject_label, "dwi", os.path.splitext(os.path.basename(dwi_file))[0] + ".bvec")
      bval = os.path.join(args.bids_dir, "sub-%s"%subject_label, "dwi", os.path.splitext(os.path.basename(dwi_file))[0] + ".bval")
      if not os.path.exists(os.path.join(args.output_dir, "dwi")):
        os.mkdir(os.path.join(args.output_dir, "dwi"))
      preproc_file = os.path.join(args.output_dir, "dwi", os.path.split(dwi_file)[-1].replace(".nii", "_preproc.mif"))
      cmd = "dwipreproc -fslgrad %s %s -rpe_none AP %s %s"%(bvec, bval, dwi_file, os.path.join(args.output_dir, preproc_file))
      print(cmd)
      # TODO remove
      cmd = "mrconvert -fslgrad %s %s %s %s"%(bvec, bval, dwi_file, preproc_file)
      subprocess.check_call(cmd, shell=True)
      if not os.path.exists(os.path.join(args.output_dir, "mask")):
        os.mkdir(os.path.join(args.output_dir, "mask"))
      mask_file = os.path.join(args.output_dir, "mask", os.path.split(dwi_file)[-1].replace("dwi.nii", "mask.mif"))
      print(cmd)
      cmd = "dwi2mask %s %s"%(preproc_file, mask_file)
      subprocess.check_call(cmd, shell=True)
      bias_file = preproc_file.replace(".mif", "_bias.mif")
      cmd = "dwibiascorrect -ants -mask %s %s %s"%(mask_file, preproc_file, bias_file)
      print (cmd)
      subprocess.check_call(cmd, shell=True)
      os.remove (preproc_file)

# running group level 1 (intensity normalisation)
elif args.analysis_level == "group1":
  mask_dir = os.path.join(args.output_dir, "mask")
  dwi_dir = os.path.join(args.output_dir, "dwi")
  norm_dwi_dir = os.path.join(args.output_dir, "norm_dwi")
  fa_template = os.path.join(args.output_dir, "fa_template.mif")
  wm_mask = os.path.join(args.output_dir, "wm_mask.mif")
  cmd = "dwiintensitynorm %s %s %s %s %s"%(dwi_dir, mask_dir, norm_dwi_dir, fa_template, wm_mask)
  print (cmd)
  subprocess.check_call(cmd, shell=True)
  response_dir = os.path.join(args.output_dir, "response")
  if not os.path.exists(response_dir):
    os.mkdir(response_dir)

# running participant level 2 (dwi2response)
elif args.analysis_level == "participant2":
  for subject_label in subjects_to_analyze:
    input_dwi = os.path.join(args.output_dir, "norm_dwi", "sub-" + subject_label + "_dwi.mif")
    output_response = os.path.join(args.output_dir, "response", subject_label + "-response.txt")
    cmd = "dwi2response tournier %s %s"%(input_dwi, output_response)
    print (cmd)
    subprocess.check_call(cmd, shell=True)

# running group level 2 (average response)
elif args.analysis_level == "group2":
  cmd = "average_response %s/* %s"%(os.path.join(args.output_dir, "response"), os.path.join(args.output_dir, "average_response.txt"))
  subprocess.check_call(cmd, shell=True)

# running participant level 3 (upsample, brain mask, CSD)
elif args.analysis_level == "participant3":
  for subject_label in subjects_to_analyze:
    input_dwi = os.path.join(args.output_dir, "norm_dwi", "sub-" + subject_label + "_dwi.mif")
    if not os.path.exists(os.path.join(args.output_dir, "upsampled_dwi")):
      os.mkdir(os.path.join(args.output_dir, "upsampled_dwi"))
    upsampled_dwi = os.path.join(args.output_dir, "upsampled_dwi", "sub-" + subject_label + "_dwi.mif")
    #TODO, detect if resolution > 2mm
    cmd = "mrresize -scale 2.0 %s %s"%(input_dwi, upsampled_dwi)
    subprocess.check_call(cmd, shell=True)
    if not os.path.exists(os.path.join(args.output_dir, "upsampled_mask")):
      os.mkdir(os.path.join(args.output_dir, "upsampled_mask"))
    input_mask = os.path.join(args.output_dir, "mask", "sub-" + subject_label + "_mask.mif")
    upsampled_mask = os.path.join(args.output_dir, "upsampled_mask", "sub-" + subject_label + "_mask.mif")
    cmd = "mrresize -scale 2.0 -interp nearest %s %s"%(input_mask, upsampled_mask)
    subprocess.check_call(cmd, shell=True)
    if not os.path.exists(os.path.join(args.output_dir, "fod")):
      os.mkdir(os.path.join(args.output_dir, "fod"))
    fod = os.path.join(args.output_dir, "fod", "sub-" + subject_label + "_fod.mif")
    group_response = os.path.join(args.output_dir, "response", subject_label + "-response.txt")
    cmd = "dwiextract %s - | dwi2fod msmt_csd - %s %s -mask %s"%(upsampled_dwi, group_response, fod, upsampled_mask)
    subprocess.check_call(cmd, shell=True)

# running group level 3 (generate study specific FOD template)
elif args.analysis_level == "group3":
  # TODO if more than 20 subjects randomly select subset to generate template
  fod = os.path.join(args.output_dir, "fod")
  upsampled_mask = os.path.join(args.output_dir, "upsampled_mask")
  fod_template = os.path.join(args.output_dir, "fod_template.mif")
  cmd = "population_template %s -mask_dir %s %s"%(fod, upsampled_mask, fod_template)
  subprocess.check_call(cmd, shell=True)

# running participant level 3 (register FODs, warp FODs, compute FD, reorient fixels, compute FC, compute FDC)
elif args.analysis_level == "participant4":
  if not os.path.exists(os.path.join(args.output_dir, "warp")):
    os.mkdir(os.path.join(args.output_dir, "warp"))
  if not os.path.exists(os.path.join(args.output_dir, "warped_fod")):
      os.mkdir(os.path.join(args.output_dir, "warped_fod"))
  fod_template = os.path.join(args.output_dir, "fod_template.mif")
  for subject_label in subjects_to_analyze:
    fod = os.path.join(args.output_dir, "fod", "sub-" + subject_label + "_fod.mif")
    upsampled_mask = os.path.join(args.output_dir, "upsampled_mask", "sub-" + subject_label + "_mask.mif")
    subject2template_warp = os.path.join(args.output_dir, "warp", "sub-" + subject_label + "subject2template_warp.mif")
    template2subject_warp = os.path.join(args.output_dir, "warp", "sub-" + subject_label + "template2subject_warp.mif")
    cmd = "mrregister %s -mask1 %s %s -nl_warp %s %s"%(fod, upsampled_mask, fod_template, subject2template_warp, template2subject_warp)
    subprocess.check_call(cmd, shell=True)
    warped_fod = os.path.join(args.output_dir, "warped_fod", "sub-" + subject_label + "_warped_fod.mif")
    cmd = "mrtransform -noreorientation %s -warp %s %s"%(fod, subject2template_warp, warped_fod)
    subprocess.check_call(cmd, shell=True)


# running group level 3 (compute fixel template mask, compute voxel mask, tractography, SIFT, fixelcfestats)
elif args.analysis_level == "group4":
  print ("group4")
