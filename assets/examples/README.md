# Example internal-standard tables

`climgrass_OE12-3_pos_equisplash_rt.csv` is an example of a lab's own standards
table with expected retention times. It is the shipped EquiSPLASH table
(positive mode) for the 10 components whose retention time in the ClimGrass
batch OE12-3 (March 2023) is known from the Skyline integration of the
internal standards in the batch's NIST SRM 1950 injections
(`ISTDSRM1950.csv`, Samrat et al. 2025 quality control). PG, CE and MG are
left out: PG was integrated in only one of three injections, CE and MG were
not integrated.

Retention times are batch- and column-specific (PC-d7 elutes at 15.6 min in
OE12-3 and 14.7–15.8 min in other reference batches), so use your own.
Pass a table like this with `--is_targets`; the QC step then accepts a
standard only within `--is_rt_window` minutes of its expected time.
