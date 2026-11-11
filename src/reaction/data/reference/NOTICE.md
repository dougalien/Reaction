# Reference data bundled with Reaction

These files are copied unchanged from their sources on 2026-11-04 so that Reaction carries the
thermodynamic data it computes with (in the hosted app and in the exe). Reaction reads them through
`reaction.app.config.external_dir()`.

| Folder | What | Source | License |
| --- | --- | --- | --- |
| `OBIGT/` | The OBIGT thermodynamic database of CHNOSZ 2.3.0: HKF parameters of aqueous species, Maier-Kelley minerals and gases | https://github.com/jedick/CHNOSZ (`inst/extdata/OBIGT`, commit 4719c5c) | GPL-3.0 (https://www.gnu.org/licenses/gpl-3.0.html) |
| `OBIGT/Berman/` | Parameters of the Berman (1988) minerals listed in `Berman_cr.csv` | https://github.com/jedick/CHNOSZ (`inst/extdata/Berman`, commit 4719c5c) | GPL-3.0 |

The EQ3/6 database itself is not bundled. Reaction reads only its ion sizes, and those numbers are in
`reaction/data/ion_sizes_eq36.tsv` (cited there). An installed `EQ36_DB/data0.cmp` in the data folder
is still used when present.

To refresh OBIGT, run `scripts/fetch_references.ps1` and copy `obigt/` to `OBIGT/`.

## `reaction/data/pitzer.dat`

PHREEQC's Pitzer database, copied unchanged from https://github.com/usgs-coupled/phreeqc3
(`database/pitzer.dat`, commit be755855, retrieved 2026-11-11, sha256
dada3271bb67d249b99c93fef870134e9cb342fc5e6aa94796fd33cacc3383ba). PHREEQC and its databases are
public domain (U.S. Geological Survey). Reaction reads its Pitzer parameters, equilibrium constants and
molar volumes (`reaction/phreeqc_pitzer.py`); the Pitzer interaction parameters are also in
`reaction/data/pitzer_phreeqc.tsv`.
