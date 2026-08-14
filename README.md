\#SEGA biosamples psychopathology - stress



\##Aim: Relate biosample data (hair and saliva) and pupillometry (baseline pupil size \[BPS] and stimulus‑evoked pupillary response \[SEPR]) to transdiagnostic psychopathology in an autism‑enriched sample.



\## Data assembly 



Combine multiple sources: a master SEGA file (demographics, biosamples), experimental task files, eye‑tracking trial‑level data, comorbidity records, and questionnaire data (CBCL, YSR, SDQ). Calculate CBCL factor scores using published loadings (McElroy 2017) and rename variables consistently.



\## Data manipulation – Derive key variables:



* Saliva cortisol change (before vs. after experiment).
* Log‑transformed hair and saliva cortisol.
* IQ estimate (mean of subtests) and BMI.
* Average BPS and SEPR across oddball conditions; change scores.
* SDQ internalizing/externalizing raw and norm‑based percentiles (self and parent).
* Global internalizing/externalizing composites combining SDQ and YSR/CBCL.
* Data quality checks – Inspect missing data, examine biosample distributions (hair vs. salivary cortisol), assess intra‑ and inter‑assay CV, and identify outliers (e.g., implausible YSR total scores).



\## Covariate screening – Test associations of age, gender, IQ, autism diagnosis, and (for biosamples) hair mass, BMI, time of day, and measurement error with psychopathology, pupillometry, and cortisol measures to guide model adjustments.



\## Sample definition – Restrict primary analyses to participants with both self‑reported psychopathology (YSR total) and hair cortisol available (complete case analysis).



\## Data Analysis



* Bivariate and multivariate modeling
* Correlate biosamples, pupillometry, and psychopathology scales.
* Fit linear models (LM) to test effects of psychopathology on cortisol (hair and saliva change) and on pupillary measures (BPS, SEPR), adjusting for relevant covariates.
* Explore combined models where pupillometry and cortisol (and their interactions) predict internalizing/externalizing scores.
* Classification – Apply logistic regression with bootstrapped AUC estimation to evaluate how well pupillometry and cortisol distinguish clinically elevated internalizing (YSR T‑score ≥ 70).
* trial‑based pupillometry – Extract individual slopes of BPS change across the visual oddball trials using linear mixed models, then relate these slopes to psychopathology.
* Condition‑based clustering – Impute missing trial‑level pupillometry data, compute Euclidean distances, and perform hierarchical clustering to identify subgroups based on task‑condition responses; compare psychopathology profiles across clusters.

