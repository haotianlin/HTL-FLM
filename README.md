# HTL-FLM

Code for Paper [On Hypothesis Transfer Learning in Functional Linear Models](https://arxiv.org/abs/2206.04277).



## Simulation

The 'Simulation' Folder contains a synthetic data demo for the proposed TL-FLR.

* './Simulation/data_gen.R' contains the R code for generating synthetic data.
* './Simulation/TL_FLR_funcs.R' contains the R code for TL-FLR/ATL-FLR implementation.
* './Simulation/TL_FLR_demo.R' contains the R code for a demo to compare target-only baseline, TL-FLR, Detect-TL, ALT-FLR (EW), and ATL-FLR.


## Real Data Application

### Application on Finantial Dataset

The stock price data are collected from [Yahoo Finance](https://finance.yahoo.com/), and we focus on stocks whose corresponding companies have a market cap over $20$ Billion with period from 4/1/2021 to 9/30/2021. We divide the sectors based on the division criteria on [Nasdaq](https://www.nasdaq.com/market-activity/stocks/screener). 

* './Application/Financial-Application/Raw-Data' contains the original stock price data for each sector.
* './Application/Financial-Application/Processed-Data' contains the processed stock price data based on the preprocess description in Appendix G.1.


### Human Activity Recognition Using Smartphones

The dataset is available at [https://archive.ics.uci.edu/dataset/240/human+activity+recognition+using+smartphones](https://archive.ics.uci.edu/dataset/240/human+activity+recognition+using+smartphones).

* './Application/HAR-Application/HAR_processed.RData' contains the processed body signal data based on the preprocess description in Appendix G.2.



## Citation
```
@misc{lin2024hypothesis,
      title={On Hypothesis Transfer Learning of Functional Linear Models}, 
      author={Haotian Lin and Matthew Reimherr},
      year={2024},
      eprint={2206.04277},
      archivePrefix={arXiv},
      primaryClass={stat.ML}
}
```