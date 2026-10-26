# Update Subworkflow WDLs

1. Inspect what subworkflows and versions are being used by the parent workflow
```
git submodule foreach 'git describe --tags'
```

2. Delete the old preprocessed WDL from /imports/, /subworkflows/ and .gitmodules if that specific subworkflow is not needed anymore (e.g. moving from bwaMem to bwamem2 would require deleting bwaMem subworkflow)
```
rm imports/pull_<subworkflow>.wdl
git rm subworkflows/<subworkflow>
```

3. Update/Add subworkflows to the workflow as git submodules 
(Reference: [this section](https://wiki.oicr.on.ca/spaces/GSI/pages/146836087/Develop+Wrapper+Workflows+Workflow+of+Workflows#DevelopWrapperWorkflows(WorkflowofWorkflows)-Workingwithgitsubmodules))
```
# For Adding a new subworkflow run
git submodule add git@github.com:oicr-gsi/<subworkflow>.git subworkflows/<subworkflow>

# For updating an existing subworkflow run
cd subworkflows/<subworkflow>
git tag --list
git checkout <subworkflow tag to be used>
cd ../..
```

4. update subworkflow "pull" imports WDL using [gsi-wdl-tools](https://github.com/oicr-gsi/gsi-wdl-tools)
```
module load gsi-wdl-tools
generate-subworkflow-import --input-wdl subworkflows/<subworkflow>/<subworkflow>.wdl --pull-all --output-wdl-path imports/pull_<subworkflow>.wdl
```

5. Remove the old imports.zip and make a new imports.zip
```
rm imports.zip
zip -X -r -9 imports.zip imports/
```
6. Edit umiCollapse.wdl to be compatible with the updated subworkflows including adding or removing any stale import "imports/pull_<subworkflow>.wdl" lines. 

7. Run tests 
```
ssh shn.hpc.oicr.on.ca

# if you don't have the workflow, make sure to clone it from shn.hpc first, as slurm nodes do not have internet access
git clone git@github.com:oicr-gsi/umiCollapse.git

salloc --account gsi --mem=4g --time=01:00:00
cd umiCollapse

# checkout the branch/tag that needs to be tested
git checkout <branch/tag to be tested>

# run the tests
/.mounts/labs/gsi/vidarr/test_vidarr_workflow.sh
```

8. Generate a new README using [gsi-wdl-tools](https://github.com/oicr-gsi/gsi-wdl-tools)
```
module load gsi-wdl-tools
generate-markdown-readme --input-wdl-path umiCollapse.wdl > README.md
```