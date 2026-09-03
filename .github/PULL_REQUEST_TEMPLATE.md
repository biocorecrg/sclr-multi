## PR checklist

- [ ] This pull request contains a clear description of the changes and their purpose.
- [ ] Tests have been added or updated where appropriate.
- [ ] The test suite passes (`nextflow run . -profile test,docker --outdir <OUTDIR>`).
- [ ] The pipeline has been checked for unexpected warnings in debug mode (`nextflow run . -profile debug,test,docker --outdir <OUTDIR>`).
- [ ] User-facing parameter changes are reflected in the relevant documentation and schema.
- [ ] `docs/usage.md` is updated where necessary.
- [ ] `docs/output.md` is updated where necessary.
- [ ] `CHANGELOG.md` is updated where appropriate.
- [ ] `README.md` is updated where appropriate.
