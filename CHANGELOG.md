# CHANGELOG



## v0.1.0 (2025-06-26)

### Feature

* feat: Add HTML reporting with interactive tables, and config-driven options

- New script of generate_report.py and template_report.html for advanced reporting
- Added config options:
  - show_count_matrix: Show deduplicated count matrix table (default: true)
  - keep_intermediate: Control retention of intermediate files (default: false)
  - html_template: Path to HTML template
- HTML report now includes:
  - Tool version at the top
  - PEAR, FASTP, pipeline, and reference stats in a summary table
  - Interactive, searchable, paginated count matrix (DataTables.js)
  - Collapsible full pipeline log
- Table and plot options are controlled via config.json
- Added --keep_intermediate flag to pipeline CLI; cleans up intermediate files unless set
- Improved parsing and deduplication of count matrix (removes Unnamed: 0, insertions, deletions, matches)
- All new features are backward compatible and configurable ([`6b5c1ab`](https://github.com/hassansaei/CapScreen/commit/6b5c1abff76275c4930fa8c98e49a696694b5266))

### Unknown

* Update README: clarify aim, usage, input/output, and add user-specific notes and license ([`ca70e9a`](https://github.com/hassansaei/CapScreen/commit/ca70e9a59444a1c0c3ce5c41d7debf3fbc544b77))

* Merge branch &#39;master&#39; of https://github.com/hassansaei/CapScreen ([`d654bbe`](https://github.com/hassansaei/CapScreen/commit/d654bbe66eb8be5262b1d05247f8cff5f8db91c1))

* Update README: clarify aim, usage, input/output, and add user-specific notes and license ([`64008cd`](https://github.com/hassansaei/CapScreen/commit/64008cd7092b7c7a76cebcabe8c9f1edee91c573))


## v0.0.1 (2025-06-24)

### Ci

* ci(release): automate package build and semantic-release with GitHub Actions ([`52c0491`](https://github.com/hassansaei/CapScreen/commit/52c0491d2c02f51eb88115d720a0fad05da4588d))

### Fix

* fix(ci): remove npm and PyPI steps for pure Python GitHub release workflow ([`6d64177`](https://github.com/hassansaei/CapScreen/commit/6d64177e179143105e57671041b2929dee5bc88e))

* fix(ci): remove npm from release workflow to support python-only project ([`e463df3`](https://github.com/hassansaei/CapScreen/commit/e463df377322c1002134ead3989cc46b6e8f3c2c))

### Unknown

* Refactor: Unify CLI, modularize count, update logging

- Rename CapScreen.py to cli.py and move to capscreen/ as the unified CLI entry point.
- Move count.py to capscreen/scripts/count.py and refactor as a module for import and CLI use.
- Update all imports and entry points to reflect new structure.
- Update setup.py to use console_scripts entry point: &#39;capscreen = capscreen.cli:main&#39;. ([`1babcd0`](https://github.com/hassansaei/CapScreen/commit/1babcd0f2f218bc02850f40bf57d51ca3df3aa85))

* Unify logging for entire CapScreen pipeline into a single log file

- Set up a single log file (sample_name.pipeline.log) at the start of the pipeline in CapScreen.py
- Pass the unified log file path to all steps, including the count step
- Modified count.py to accept an optional log_file argument for logging
- Ensured no duplicate log handlers are added in count.py ([`2e95468`](https://github.com/hassansaei/CapScreen/commit/2e95468e4edecd18af50e6964697a015aee5b6e3))

* Ensure unassigned reads are labeled as &#39;Unassigned&#39; in ID_WLG column and count table include informative columns ([`6fdbb3d`](https://github.com/hassansaei/CapScreen/commit/6fdbb3dfec33d53bebe7df6b2dbfecb504f7aa83))

* Add pandas and biopython packages to the conda env ([`e3a53b0`](https://github.com/hassansaei/CapScreen/commit/e3a53b0ef02df830a82fa3000adb472121cad943))

* Refactor CapScreen CLI to support subcommands and modular pipeline

- Refactored CapScreen.py to use argparse subparsers for three subcommands: `pipeline`, `align`, and `count`.
- Integrated count.py as a module, allowing the `count` subcommand to run the variant counting step independently.
- The `pipeline` subcommand now runs QC, alignment, and counting in sequence.
- The `align` subcommand runs only QC and alignment.
- Logging is consistent and clear for each subcommand, with sample-specific log files for QC/alignment and separate logs for counting.
- Updated argument parsing for clarity and modularity.
- Verified that config.json is complete and consistent with the new pipeline structure. ([`fd4c0fd`](https://github.com/hassansaei/CapScreen/commit/fd4c0fda1d565b02004df59cf2c11d3b4d0905a9))

* New module: count.py to count and normalize capsid variants ([`71c335b`](https://github.com/hassansaei/CapScreen/commit/71c335bca7d2a81b99c216749a36abfd81ede3ce))

* Initial commit ([`1c9ab62`](https://github.com/hassansaei/CapScreen/commit/1c9ab62a8852c6d163e43c27e36ee9c58662655b))
