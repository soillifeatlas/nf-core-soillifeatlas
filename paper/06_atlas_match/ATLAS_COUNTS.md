# Which Atlas count belongs where (review item 5)

2026-11-25. Sources: `~/Desktop/SoilLifeAtlas/release/release.json` (release 0.9 draft, SLA `fd7e820`,
2026-10-19) and the release reprocessing on LiSC (`sla_reprocess/t3/analysis_pos/a_presence_long.csv`)
used for the soil-to-Atlas matching.

| Count | POS | NEG | What it is | Where to use it |
|---|---|---|---|---|
| Reference samples collected | 168 | 195 | every reference injection in the release | Fig. 1A, abstract, graphical abstract, §3.1 "collection" |
| Phyla collected | 19 | 19 | includes 3 phyla with a single sample each: Cyanobacteriota, Mortierellomycota, Cercozoa | Fig. 1A, §3.1 |
| Analysis samples | 164 | 192 | samples kept for analysis (release `analysis_samples_*`) | §2.1 / §3.1 when describing the analysed reference |
| Phyla analysed | 16 | 16 | phyla with replication (rule: n ≥ 2 in both polarities); the 3 single-sample phyla are dropped | Fig. 1A note, §3.1, Atlas occurrence / enrichment statements |
| Reference samples with ≥ 1 detected feature in the reprocessing | 161 | not used | denominators for per-phylum detection in matching | Fig. 6C, §3.6 matching ("161 reference samples, 16 phyla") |
| Group sizes (POS) | Bacteria 50, Fungi 45, Archaea 19, Plants 19, Protists 18, Soil animals 17 | Bacteria 56, Fungi 51, Plants 26, Soil animals 27, Archaea 21, Protists 14 | samples collected per group | graphical abstract |
| Atlas lipids | 11,371 | 5,697 | release lipid counts | §3.1, Fig. 1B caption |

161 vs 164 (POS): the reprocessing presence table has 3 fewer samples than the release's analysis
set (per phylum: Mollusca 1 vs 2, Nematoda 3 vs 4, Thermoproteota 2 vs 3). These are the samples
that contribute no detected feature to the presence table; the matching therefore uses 161.

Fig. 1A ("168 samples, 19 phyla") is correct for the **collected** POS reference. Statements about
occurrence, matching or enrichment should use **16 analysed phyla** (and 164 POS analysis samples,
or 161 for the matching denominators). All release 0.9 numbers are still marked draft; freeze them
at release 1.0 (issue #12).
