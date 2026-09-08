# Excluded code policy

Historical duplicates and exploratory scripts remain in the original workspace but are intentionally absent here. Exclusion does not imply that the original files are invalid; it means that no direct link to a reported manuscript panel was established during this audit.

Excluded categories include:

- v11 and v12 benchmark directories superseded by the selected v13 master;
- files with `old`, `old1`, `old2`, `old3`, `副本`, or `已恢复` in their names;
- unused newplot effect-size and three-set forest outputs;
- cached numerical results, iteration checkpoints, plot queues, and temporary reports;
- generated package documentation, websites, IDE state, and R CMD check directories;
- upstream mapping scripts, because the manuscript Methods provides the processing software and parameters and the author chose not to distribute those local workflow wrappers;
- raw or processed study data, which should be deposited separately and referenced by accession or stable archive identifier.

If a later figure audit establishes that an excluded script produced a submitted panel, copy that exact source file into this snapshot and add one row to `source_provenance.tsv` and `figure_script_map.tsv`. Do not silently replace an existing copied script.
