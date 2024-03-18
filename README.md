# A python API to retrieve MOF information from FAIR-MOFs database
This API enables users to access all the data we uploaded into NOMAD. The data set is a subset of 35k thoroughly curated MOFs structures that maps each MOF to its building units and porosity property.  Notably, the data comes alongside tools that we implemented in NOMAD to facilitate the search of MOFs from their building units and porous property.

The software (porosityNormaliser) is available to all users to upload cif files from experiment and output from various computational chemistry packages. The software will identify the structure as either MOF, COF, zeolite or a porous system if the pore limiting diameter of the system is greater than 1.8 Å.

On the other hand, the MOF dataset that we uploaded contains geometry optimised structures of curated MOFs using the GFN-xTB method. A second upload containing  experimental synthetic conditions extracted from journal articles will be published in December.

This is our endeavour to create a comprehensive and FAIR  MOF  database that maps the structure of all experimentally synthesised MOFs to their building units, properties and experimental synthetic conditions. Our overarching goal is to design new machine learning models for predicting the synthesis conditions of any already synthesised or hypothetical MOF.

To enable easy access to this database, this api will enable users to directly extract the data from NOMAD without needing to worry about the REST API. Here all you need is the CSD refcode for the MOF, which is not available for the data set in NOMAD.

# Installation
## Step 1:clone the repository
Start by cloning the repository the git repository
``` git clone https://github.com/bafgreat/fairmof-api.git ```
## Step 2: install
Go to the respository and pip install by following commands below
i. `cd fairmof-api`
ii. `pip install .`

At this stage the api will be installed and you can now use it as a library to download any mof using the csd refcode.

# Usage
This module is still in development and more functions will be added in the future.
## api library
At the moment you can download any mof, geometric properties and building units using the following commands.

### Import the package
Import the package with the following command
``` from fairmof_api.api import download_archive ```

### Download the mof using the following command
 You can download any MOF whose recode is present using the following commadn download_archive(refcode), where refcode refers to the csd identifier. E.g.
 ```download_archive("EDUSIF")```
 The above command will download MOF-5 whose csd refcode is EDUSIF. This will create a folder called FAIR-MOFs in the same working directory. In this directory, 2 json files will be create and two folders will be created.

#### mof_properties.json
 The file called "mof_properties.json" contains the porosity properties.
#### mof_sbu.json
The file "mof_sbu.json" contains all information about the building units.

#### mofs_recognised_by_nomad
The folder "mofs_recognised_by_nomad" contains the cif file of the MOF if NOMAD recognizes this system as a MOF. It is important to note that NOMAD recognizes any system as a MOF if it is periodic, have atleas one metal and organic atoms with  a pore limiting diameter greater than 1.8 Å.

#### mofs_not_recognised_by_nomad
The folder "mofs_not_recognised_by_nomad" contains cif file of the refcode in case csd recognises it as a MOF but nomad does not.

# N.B.
You can also perform a batch operation, where every information will still be appended in FAIR-MOFs folder. You could aslo provide a custom name for the result folder.

## Running batch operations
``` from fairmof_api.api import download_archive ```
```
for refcode in list_of_refcodes:
    download_archive(refcode)
```
The output of the above operation will be written to in FAIR-MOFs folder. if you do not want it to be written to this folder you can provide a custom name or path to the folder.
``` download_archive(refcode, result_folder="Path_to_folder", extension='cif') ```
Notice that you can also chose the extension. if you do not want the output to be written to be in cif format. You can choose any ase writable extension. 








