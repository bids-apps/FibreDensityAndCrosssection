### Description
This BIDS App enables group analysis of diffusion MRI data by performing a Fixel-Based Analysis (FBA) of Fibre Density, Fibre Cross-section and a combined measure (Fibre Density & Cross-section).

The analysis pipeline relies primarily on the [MRtrix3](www.mrtrix.org) software package.

### Documentation

Full documentation for MRtrix3 is available online [here](http://userdocs.mrtrix.org/). This pipeline performs pre-processing of diffusion MRI data according to the steps outlined [here](http://mrtrix.readthedocs.io/en/latest/workflows/DWI_preprocessing_for_quantitative_analysis.html), and fixel-based analysis using the steps listed [here](http://mrtrix.readthedocs.io/en/latest/workflows/fixel_based_analysis.html).

### Error Reporting

For help and support please post a question on the [MRtrix3 discussion](http://community.mrtrix.org/).

### Acknowledgement

When using this pipeline, please include the following paragraph to descibe the method used. Citations can be found in the attached [bibtex](./fixel-based_analysis.bib) file.

Pre-processing was performed by first denoising the diffusion-weighted images (Veraart et al. 2016), followed by eddy-current and motion correction (Andersson et al. 2016). Fibre orientation distributions (FOD) were computed using multi-tissue constrained spherical deconvolution (Jeurissen et al. 2014). Simultaneous bias field correction and intensity normalisation (across subjects) was performed as per Raffelt et al. 2017. All subjects were registered to a study-specific FOD template using an iterative update approach (Raffelt et al. 2011, Raffelt et al. 2012b). Three quantitative measures were computed for each white matter fixel: apparent Fibre Density (FD) (Raffelt et al. 2012), Fibre Cross-section (FC)(Raffelt et al. 2017) and also a combined measure Fibre Density and Cross-section (FDC)(Raffelt et al. 2017b).Statistical analysis was performed using connectivity-based fixel enhancement (CFE)(Raffelt et al. 2015). CFE exploits fixel-fixel connectivity information derived from whole-brain fibre tractography streamlines computed on the FOD template (Tournier et al. 2010). Following tractography, SIFT was used to reduce tractogram biases (Smith et al. 2013). We assigned family-wise error corrected p-values to each fixel using permutation testing of the CFE enhanced t-statistics (5000 permutations).

If your BIDS dataset has reverse phase encoded b=0 or DWI image pairs, please also add this sentence to the paragraph above (after the eddy-current correction): Suceptibility-induced distortions were corrected with reverse phase encoded pairs (Andersson et al. 2003).


### Instructions

This pipeline requires that data be organized in accordance with the [BIDS](http://bids.neuroimaging.io) spec.


**To get your container ready to run type the following your terminal:**
```{bash}
$ docker pull bids/fibredensityandcrosssection
```

**Before starting, let's check out the help page**

```
$ docker run -ti bids/fibredensityandcrosssection -h

SYNOPSIS

     run.py [ options ] bids_dir output_dir analysis_level

        bids_dir     The directory with the input dataset formatted according to
                     the BIDS standard.

        output_dir   The directory where the output files should be stored. If
                     you are running group level analysis this folder should be
                     prepopulated with the results of the participant level
                     analysis.

        analysis_level Level of the analysis that will be performed. Valid
                     choices are: [participant1, group1, participant2, group2,
                     participant3, group3, participant4, group4].  Multiple
                     participant level analyses can be run independently(in
                     parallel) using the same output_dir.

DESCRIPTION

     Perform group analysis of diffusion MRI data with a Fixel-Based Analysis
     (FBA) of Fibre Density, Fibre Cross-section and a combined measure (Fibre
     Density & Cross-section). The analysis pipeline relies primarily on the
     MRtrix3 software package (www.mrtrix.org).

Options for the population_template script

  --participant_label PARTICIPANT_LABEL
     The label(s) of the participant(s) that should be analyzed. The label
     corresponds to sub-<participant_label> from the BIDS spec (so it does not
     include "sub-"). If this parameter is not provided all subjects should be
     analyzed. Multiple participants can be specified with a space separated
     list.

  --n_cpus INT
     The number of CPU cores available on the compute node. Set to 0 to use the
     maximum number of cores available

Standard options

  -continue <TempDir> <LastFile>
     Continue the script from a previous execution; must provide the temporary
     directory path, and the name of the last successfully-generated file

  -force
     Force overwrite of output files if pre-existing

  -help
     Display help information for the script

  -nocleanup
     Do not delete temporary files during script, or temporary directory at
     script completion

  -tempdir /path/to/tmp/
     Manually specify the path in which to generate the temporary directory

  -quiet
     Suppress all console output during script execution

  -verbose
     Display additional information for every command invoked

AUTHOR
     David Raffelt (david.raffelt@florey.edu.au)

```

**To run the pipeline using docker**

This fixel-based analysis pipeline has been broken up into several stages, each defined by the "analysis level" positional argument. Each level is labelled as either participant or group. Participant levels can be run on different subjects independently, while group level analysis is performed on all subjects within the group. The order in which the analysis should be run is participant1, group1, participant2, group2, participant3, group3, particpant4, group4.

In order to share data between our container and the rest of our machine, we need to mount a volume. Docker does this with the `-v` flag. Docker expects its input formatted as: `-v path/to/local/data:/path/in/container`.

For example, to run the first particpant level analysis on a single subject use:

```{bash}

docker run -i --rm -v /path/to/local/data:/bids_input -v /path/to/local/output:/output bids/fibredensityandcrosssection /bids_input /output participant1 --participant_label 01
```

For example, to run the first particpant level analysis on a single subject use:
