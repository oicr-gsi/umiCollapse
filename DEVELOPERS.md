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

4. Replace the old WDL file in /subworkflows/ with the new one
```
cp path/to/subworkflow.wdl subworkflows/subworkflow.wdl
```

5. Preprocess the new WDL with dockstore_preprocess.py found in [gsi-wdl-tools](https://github.com/oicr-gsi/gsi-wdl-tools)
```
python3 path/to/gsi-wdl-tools/dockstore_preprocess.py --input-wdl-path path/to/umiCollapse/subworkflows/subworkflow.wdl --pull-all True --output-wdl-path path/to/umiCollapse/imports/pull_subworkflow.wdl
```

6. Make a new imports.zip
```
zip -r path/to/umiCollapse/imports.zip path/to/umiCollapse/imports/
```

7. Submit to Cromwell for testing
```
java -jar $cromwell submit umiCollapse.wdl --inputs path/to/test.json --imports imports.zip --host http://cromwell-dev.hpc.oicr.on.ca:8000
```

8. Generate a new README using generate-markdown-readme.py also found in [gsi-wdl-tools](https://github.com/oicr-gsi/gsi-wdl-tools)
```
python3 path/to/gsi-wdl-tools/generate_markdown_readme.py --input-wdl-path path/to/umiCollapse/umiCollapse.wdl > README.md
```