```diff
---------  Under construction!  ---------
```
# mutationsignatures-nf
De novo mutational signature as described by Alexandrov et al., 2015
<br><br>

## Quick Start

Required Arguments:

| argument       | value | 
|:--------------:|:-----:| 
| `maf_folder`| a path to the input folder with `.maf` files from individuals to be merged to create the cohort `.maf`| 

To test the pipeline with the example input you can run:

```nextflow
# Clone the repository
git clone https://github.com/cgpu/mutationalsignatures-nf.git

# cd into the repo folder 
cd mutationalsignatures-nf/

# Execute nextflow run command with example input parameters
nextflow run mutationalsignatures-nf/main.nf --maf_folder path/to/input/maf_folder/ 
```

