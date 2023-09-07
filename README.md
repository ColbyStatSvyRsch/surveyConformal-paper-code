# surveyConformal-paper-code

R code and output associated with "Design-based conformal prediction," Wieczorek (2023+), accepted to [*Survey Methodology*](https://www150.statcan.gc.ca/n1/en/catalogue/12-001-X) for publication in late 2024. Preprint: https://arxiv.org/abs/2303.01422



## Data

The examples in Section 4.1 rely on data from the Medical Expenditure Panel Survey (MEPS). We do not share the raw data here. It is available for download from the Agency for Healthcare Research and Quality (AHRQ):  
https://meps.ahrq.gov/mepsweb/index.jsp

MEPS data users must follow the Data Use Agreement:  
https://meps.ahrq.gov/data_stats/data_use.jsp


The simulations in Section 4.2 rely on data from the Academic Performance Index (API). We accessed that data through [the `{survey}` R package](https://cran.r-project.org/package=survey), version 4.1-1.
