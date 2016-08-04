This BIDS App enables group analysis of diffusion MRI data by performing a Fixel-Based Analysis (FBA) of Fibre Density, Fibre Cross-section and a combined measure (Fibre Density & Cross-section).

The analysis pipeline relies primarily on the [MRtrix3](www.mrtrix.org) software package.

### Documentation

Full documentation for MRtrix3 is available online [here](http://userdocs.mrtrix.org/). This pipeline performs pre-processing of diffusion MRI data according to the steps outlined [here](http://mrtrix.readthedocs.io/en/latest/workflows/DWI_preprocessing_for_quantitative_analysis.html). An overview of the steps performed in this fixel-based analysis are described [here](http://mrtrix.readthedocs.io/en/latest/workflows/fixel_based_analysis.html).

### Error Reporting

For help and support please post a question on the MRtrix3 discussion forum found [here](http://community.mrtrix.org/).

### Acknowledgement

When using this pipeline, please include the following paragraph to descibe the method used. Citations can be found in the attached [bibtex](./fixel-based_analysis.bib) file.


### Instructions

This pipeline requires that data be organized in accordance with the [BIDS](http://bids.neuroimaging.io) spec.


**To get your container ready to run just follow these steps:**

- In your terminal, type:
```{bash}
$ docker pull bids/FixelAnalysis
```

**Before starting, let's check out the help page**

```{bash}
$ docker run -ti bids/FixelAnalysis -h

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
This fixel-based analysis pipeline has been broken up into several stages, each defined by the "analysis level" positional argument. Each level is labelled as either participant or group. Participant levels can be run on different subjects independently, while group level analysis is performed on all subjects within the group. The order in which the analysis should be run is participant1, group1, participant2, group2, participant3, group3, particpant4, group4.

In order to share data between our container and the rest of our machine, we need to mount a volume. Docker does this with the `-v` flag. Docker expects its input formatted as: `-v path/to/local/data:/path/in/container`.

For example, to run the first particpant level analysis on a single subject use:

```{bash}

docker run -i --rm -v /path/to/local/data:/bids_input -v /path/to/local/output:/output bids/FixelAnalysis /bids_input /output participant1 --participant_label 01
```


