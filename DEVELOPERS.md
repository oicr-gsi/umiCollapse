# Update Subworkflow WDLs

1. Working with the repository (as it has git submodules)
```
# cloning
git clone --recurse-submodules <url>

# pulling changes
git pull --recurse-submodules
```

2. Inspect what subworkflow versions are used 
```
$ git submodule foreach 'git describe --tags'
Entering 'subworkflows/bamQC'
5.2.0
Entering 'subworkflows/bwamem2'
1.1.2
```

3. Updating subworkflow versions
(Reference: [Working with git submodules](https://wiki.oicr.on.ca/spaces/GSI/pages/146836087/Develop+Wrapper+Workflows+Workflow+of+Workflows#DevelopWrapperWorkflows(WorkflowofWorkflows)-Workingwithgitsubmodules))
```
# checkout the version you want to update, for example:
cd subworkflows/bamQC
git checkout 5.x.x
```

4. Building imports.zip
```
# then update "pull" and imports.zip to update what will be run
module load gsi-wdl-tools
rm -rf imports imports.zip
mkdir -p imports
for wdl in subworkflows/*/*.wdl; do
  dir=$(dirname "$wdl")
  echo "$dir: $(git -C "$dir" describe --tags)"
  generate-subworkflow-import \
    --input-wdl "$wdl" \
    --pull-all \
    --output-wdl-path "imports/pull_$(basename "$wdl")"
done
zip -X -r -9 imports.zip imports/
```

5. Edit umiCollapse.wdl to be compatible with the updated subworkflows including adding or removing any stale import "imports/pull_<subworkflow>.wdl" lines. 

6. Run tests 
```
ssh shn.hpc.oicr.on.ca

# Note: if you don't have the workflow, make sure to clone it from shn.hpc first, as slurm nodes do not have internet access

salloc --account gsi --mem=4g --time=01:00:00
cd umiCollapse

# checkout the branch/tag that needs to be tested
git checkout <branch/tag to be tested>

# run the tests
/.mounts/labs/gsi/vidarr/test_vidarr_workflow.sh
```

7. Generate a new README using [gsi-wdl-tools](https://github.com/oicr-gsi/gsi-wdl-tools)
```
module load gsi-wdl-tools
generate-markdown-readme --input-wdl-path umiCollapse.wdl > README.md
```