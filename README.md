# surveyConformal-paper-code

R code and output associated with ["Design-based conformal prediction,"](https://www150.statcan.gc.ca/n1/pub/12-001-x/2023002/article/00007-eng.htm) (Wieczorek, *Survey Methodology*, 49 (2), 2023).  
Preprint version: https://arxiv.org/abs/2303.01422



## Data

The examples in Section 4.1 rely on data from the Medical Expenditure Panel Survey (MEPS). We do not share the raw data here. It is available for download from the Agency for Healthcare Research and Quality (AHRQ):  
https://meps.ahrq.gov/mepsweb/index.jsp

MEPS data users must follow the Data Use Agreement:  
https://meps.ahrq.gov/data_stats/data_use.jsp


The simulations in Section 4.2 rely on data from the Academic Performance Index (API). We accessed that data through [the `{survey}` R package](https://cran.r-project.org/package=survey), version 4.1-1.
