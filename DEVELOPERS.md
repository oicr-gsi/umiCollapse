# Update Subworkflow WDLs

1. Edit umiCollapse.wdl to be compatible with the updated subworkflows

2. Delete the old imports.zip
```
rm imports.zip
```

3. Delete the old preprocessed WDL from /imports/ (prefixed with "pull_")
```
rm imports/pull_subworkflow.wdl
```

4. Update the subworkflow submodule to the new version 
(see [this section](https://wiki.oicr.on.ca/spaces/GSI/pages/146836087/Develop+Wrapper+Workflows+Workflow+of+Workflows#DevelopWrapperWorkflows(WorkflowofWorkflows)-Workingwithgitsubmodules) for more details on updating an existing subworkflow)
```
git submodule add git@github.com:oicr-gsi/subworkflow.git subworkflows/subworkflow
```

5. update subworkflow "pull" imports WDL using [gsi-wdl-tools](https://github.com/oicr-gsi/gsi-wdl-tools)
```
generate-subworkflow-import --input-wdl subworkflows/subworkflow/subworkflow.wdl --pull-all --output-wdl-path imports/pull_subworkflow.wdl
```

6. Make a new imports.zip
```
zip -r -9 path/to/umiCollapse/imports.zip path/to/umiCollapse/imports/
```

7. Submit to Cromwell for testing
```
java -jar $cromwell submit umiCollapse.wdl --inputs path/to/test.json --imports imports.zip --host http://cromwell-dev.hpc.oicr.on.ca:8000
```

8. Generate a new README using generate-markdown-readme.py also found in [gsi-wdl-tools](https://github.com/oicr-gsi/gsi-wdl-tools)
```
python3 path/to/gsi-wdl-tools/generate_markdown_readme.py --input-wdl-path path/to/umiCollapse/umiCollapse.wdl > README.md
```