# ORT: Omics to Reactive Transport
This repository complements the manuscript [ORT: A workflow linking genome-scale metabolic models with reactive transport codes](https://www.biorxiv.org/content/10.1101/2021.03.02.433463v1) and the KBase narrative [here](https://narrative.kbase.us/narrative/71260) which describes a workflow for linking omics data to reactive transport simulations.

Files included here:
* `ort.py` - script for downloading and translating KBase results to a plaintext file that can be used to automatically (or manually) generate PFLOTRAN infiles
* `config.yml` - configuration file with information needed to connect to the narrative associated with this work. **Note: you will need to substitute in your actual auth token and narrative details to use this with another narrative**.
* `nitrogen_cycling_0d.in` - PFLOTRAN infile for literature-based model
* `nitrogen_cycling_in_hanford_300_river_sediments_0d.in` - PFLOTRAN infile for KBase-derived model

## Configuration File
The configuration yaml contains the following:
```
endpoint: https://kbase.us/services
workspace-url: https://kbase.us/services/ws
KB_AUTH_TOKEN: <YOUR_AUTH_TOKEN_HERE>
NAR_ID: rlrubinstein:narrative_1599076553976
NAR_NAME: ORT_N_cycle
```

`endpoint` and `workspace-url` can likely be used as they are, but if you have trouble you should check the [KBase documentation](https://kbase.github.io/kb_sdk_docs/howtos/workspace.html) for updates.

`KB_AUTH_TOKEN` should be the token for an account with access to the narrative you want to access.

`NAR_ID` is the narrative ID for the narrative you want to access. To use this script with a different narrative, substitute the appropriate narrative ID.

`NAR_NAME` is the display name you want to use for storing the results. You can put whatever you want here.

To use this with narratives other than the test case shown here,

## Using this code
To use this code, you need to populate the `config.yml` file with the appropriate information. Then, run 

```
python ort.py
```

from the terminal and it will automatically download and format the reactions. A new folder will be created inside the `pulled` directory for each narrative. The narrative-specific folder will contain a `*.txt` file with a reaction string and yield term for each model in the narrative. The `*.dat` files may be used to augment the PFLOTRAN geochemical database file if the one you have does not already contain the necessary compounds. Because we left the unedited version of the Nitrososphaeraceae model in the narrative, we end up with four reactions in the txt file and four `*.dat` files even though we only use three in our PFLOTRAN model. Note that the geochemical database included with the infiles will work for these infiles, but if you use this script with your own narrative, you may end up with some compounds not currently included.

The PFLOTRAN infiles included in this repository (in the `infiles` folder) are for the 0D system described in the manuscript (generated programmatically), but the same process could be used to generate MICROBIAL\_REACTION inputs for any model you wish.

If you do not have PFLOTRAN installed, you can run these infiles using [our free PFLOTRAN v3 docker image](https://hub.docker.com/r/subsurfaceinsights/pflotran)

