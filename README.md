# Propensity score weighting plus regression adjustment does not equal doubly robust
## Reproducibility materials

### What is this?

This contains a series of R scripts to reproduce the simulation results in the above-titled paper. 

### How do I run this code? 

#### Prerequisites
You will need a recent version of R, along with the following packages `c("stdReg", "ggplot2", "data.table", "parallel", "sandwich")`. 

This uses the `targets` package for reproducibility. Read more about how to use it here: https://books.ropensci.org/targets/

#### Steps

1. Download the repository
2. Source the `_targets.R` script. Adjust the parameters `B` and `cores` in the 
3. Run `tar_make()` or `tar_make_clustermq()` to run the simulation. Warning, it takes a long time because each simulation replicate includes bootstrapping. 
4. Inspect the results with `tar_read(linear_summary)`. Also there is a pdf generated for the example analysis. 

```{r}
source("_targets.R")
tar_make()

tar_read(sim_summary)
tar_read(funk_compare)
tar_read(odd_text)
```


### Visualization

```mermaid
graph LR
  subgraph legend
    direction LR
    x0a52b03877696646([""Outdated""]):::outdated --- x7420bd9270f8d27d([""Up to date""]):::uptodate
    x7420bd9270f8d27d([""Up to date""]):::uptodate --- xbf4603d6c2c2ad6b([""Stem""]):::none
  end
  subgraph Graph
    direction LR
    xda319e07d75bcf17(["simulation_results"]):::outdated --> x389ac23a8bc5790a(["sdev_table"]):::outdated
    xda319e07d75bcf17(["simulation_results"]):::outdated --> x528a48b3a0beaa4a(["funk_compare"]):::outdated
    xda319e07d75bcf17(["simulation_results"]):::outdated --> x096535d843ac58b4(["odd_text"]):::outdated
    xda319e07d75bcf17(["simulation_results"]):::outdated --> x8df1cfa1cec0f4ca(["sim_summary"]):::outdated
    x9e69104330fdf6f8(["simulate_inverse_gaussian_inverse_gaussian_weighted_standardized_2000_200"]):::uptodate --> xda319e07d75bcf17(["simulation_results"]):::outdated
    x71467bd26b2548ad(["simulate_linear_linear_compare_100_2"]):::outdated --> xda319e07d75bcf17(["simulation_results"]):::outdated
    xe01042d5af18e30e(["simulate_linear_linear_compare_1000_2"]):::outdated --> xda319e07d75bcf17(["simulation_results"]):::outdated
    x116f146dea47cfe7(["simulate_linear_linear_compare_2000_2"]):::outdated --> xda319e07d75bcf17(["simulation_results"]):::outdated
    x9eae9ce78ce473f4(["simulate_linear_linear_compare_500_2"]):::outdated --> xda319e07d75bcf17(["simulation_results"]):::outdated
    x00130b7e30e0600b(["simulate_linear_ols_weighted_2000_2"]):::uptodate --> xda319e07d75bcf17(["simulation_results"]):::outdated
    x9d6fa315a53f6e26(["simulate_linear_ols_weighted_standardized_2000_2"]):::uptodate --> xda319e07d75bcf17(["simulation_results"]):::outdated
    x1e015faaad899be1(["simulate_linearodd_ols_weighted_standardized_odd_2000_2"]):::outdated --> xda319e07d75bcf17(["simulation_results"]):::outdated
    x5ea0a7387f8b2ef9(["simulate_log_poisson_poisson_weighted_standardized_2000_2"]):::uptodate --> xda319e07d75bcf17(["simulation_results"]):::outdated
    xac3d6e69f58e5f12(["simulate_logit_binomial_log_binomial_weighted_standardized_2000_20"]):::uptodate --> xda319e07d75bcf17(["simulation_results"]):::outdated
    x99723e3e80a40ee6(["simulate_logit_binomial_logit_binomial_weighted_standardized_2000_2"]):::uptodate --> xda319e07d75bcf17(["simulation_results"]):::outdated
    x4cc6fc8d5fa9053a(["simulate_logit_binomial_logit_compare_100_5"]):::outdated --> xda319e07d75bcf17(["simulation_results"]):::outdated
    x7649e25e048e52a8(["simulate_logit_binomial_logit_compare_1000_5"]):::outdated --> xda319e07d75bcf17(["simulation_results"]):::outdated
    xee6f069e485dd697(["simulate_logit_binomial_logit_compare_2000_5"]):::outdated --> xda319e07d75bcf17(["simulation_results"]):::outdated
    xe7cf092df7d54f75(["simulate_logit_binomial_logit_compare_500_5"]):::outdated --> xda319e07d75bcf17(["simulation_results"]):::outdated
    xe6eda53558c41c5e(["example"]):::outdated --> xe6eda53558c41c5e(["example"]):::outdated
  end
  classDef outdated stroke:#000000,color:#000000,fill:#78B7C5;
  classDef uptodate stroke:#000000,color:#ffffff,fill:#354823;
  classDef none stroke:#000000,color:#000000,fill:#94a4ac;
  linkStyle 0 stroke-width:0px;
  linkStyle 1 stroke-width:0px;
  linkStyle 21 stroke-width:0px;
```
