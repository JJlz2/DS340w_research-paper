# Data

| File | What it is | Size | Source |
|---|---|---|---|
| `raw/AK_ECC_datasheet_Kojic.xlsx` | **The parent paper's dataset.** One row per outcome measure per study: group pre/post means and SDs, change scores, sample sizes, pre-post correlations, and study characteristics (training status, duration, volume matching, exercise type, risk of bias). | 27 rows × 48 columns, 9 studies | Amdi & King (2025), OSF: https://osf.io/s6uqn/ (identical to `s6uqn-osfstorage-archive/AK ECC datasheet Kojic.xlsx`) |
| `raw/Hao2026_pone.0354725.s001.xlsx` | **Raw participant data from Hao et al. (2026)**, a trial published after the parent paper's search, used to test the parent paper's predictions. Pre/post squat 1RM, CMJ, squat jump, 30 m sprint, body mass and height for every participant. Also holds a sheet of training repetitions and the authors' SPSS output. | 30 participants, 10 per group (2 s, 4 s and self-selected eccentric tempo) | PLOS One supplementary file S1, CC BY 4.0, doi:10.1371/journal.pone.0354725 |

## Train / test / validation splits
`splits/` divides the parent paper's datasheet into three sets. The columns are the same as in `raw/AK_ECC_datasheet_Kojic.xlsx`. Whole studies are assigned to a set, never single rows, because rows from the same study share participants: splitting a study across sets would leak it into its own test.

| File | Studies | Rows (strength / hypertrophy / CMJ) |
|---|---|---|
| `splits/train_dataset.xlsx` | Dias et al. (2015), Pereira et al. (2016), Kojic et al. (2021), Pearson et al. (2022), Azevedo et al. (2022), Segers et al. (2022) | 15 (7 / 7 / 1) |
| `splits/test_dataset.xlsx` | Mike et al. (2017), Kojic et al. (2025) | 7 (2 / 4 / 1) |
| `splits/validation_dataset.xlsx` | Shibata et al. (2021), the only study with strength, hypertrophy and CMJ data | 5 (1 / 3 / 1) |

That is 6 / 2 / 1 of the 9 studies (67% / 22% / 11%).

The parent paper's Tables 1–3 are also saved as CSV in `../figures and data/`.
