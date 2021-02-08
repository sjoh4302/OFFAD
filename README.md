# OffPeriodDetection

# What is EEGLAB?
EEGLAB is an open source signal processing environment for electrophysiological signals running on Matlab and Octave (command line only for Octave). This folder contains original Matlab functions from the EEGLAB (formerly ICA/EEG) Matlab toolbox, all released under the Gnu public license (see eeglablicence.txt). See the EEGLAB tutorial and reference paper (URLs given below) for more information.

# Installing/cloning
**Recommended:** Download the official EEGLAB release from https://sccn.ucsd.edu/eeglab/download.php

**Do not download a ZIP file directly from GIT as it will not contain EEGLAB submodules**. Instead clone the reposity while pulling EEGLAB sub-modules.

```
git clone --recurse-submodules https://github.com/sccn/eeglab.git
```

If you forgot to clone the submodule, go to the eeglab folder and type

```
git submodule update --init --recursive --remote
git pull --recurse-submodules
```
