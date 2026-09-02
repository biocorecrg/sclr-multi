# Contributing to scLR-multi

Contributions to scLR-multi are welcome, including bug fixes, improvements to existing workflows, support for additional single-cell long-read technologies, documentation improvements, and new or updated tests.

## Issues

Before making substantial changes, please open an issue describing the proposed change. Bug reports should include sufficient information to reproduce the problem, including the scLR-multi version, Nextflow version, command and parameters used, and relevant error messages.

## Development

To contribute code:

1. Fork the repository and create a branch for your changes.
2. Follow the existing Nextflow DSL2 structure and coding conventions.
3. Reuse modules from [nf-core/modules](https://github.com/nf-core/modules) where appropriate rather than duplicating existing functionality.
4. Add or update tests when modifying pipeline functionality.
5. Update the relevant documentation and parameter definitions when adding or changing user-facing functionality.
6. Submit a pull request describing the changes and their purpose.

Please keep pull requests focused on a single feature or issue where possible.

## Testing

Changes affecting pipeline functionality should be tested using the appropriate test inputs and configurations described in the scLR-multi README and provided in the `test_input/` directory.

## Acknowledgement of nf-core

scLR-multi was initially derived from [nf-core/scnanoseq](https://github.com/nf-core/scnanoseq) and continues to use components and infrastructure developed by the [nf-core](https://nf-co.re) community.
