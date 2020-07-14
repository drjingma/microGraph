# microGraph
Graphical models for microbiome data

 * Fix number of variables $p=200$ and vary the sample size $n=100,200,500$.
 
 * Run both null model (a. shuffle an available dataset, e.g. amgut1.filt from *SpiecEasi* package; b. generate from Dirichlet Multinomial) and parametric/alternative model (a. a log-normal distribution; b. generative copula model used in SPIEC-EASI). 
 
      - To simulate from Dirichlet Multinomial, we can use 
 
      - For log-normal distribution, generate network (covariance or inverse covariance) from a random graph. You can use the ones I used in the code example.
      
      - For copula model, we can follow examples in SPIEC-EASI.
 
 * Implement methods including: Conet, SparCC, CCLasso, COAT, gCoda, SPIEC-EASI, SPRING. Note SPIEC-EASI is essentially running graphical lasso on CLR transformed data. Given a CLR transformed data, we can also apply other methods developed for Gaussian data, such as [SpaceJAM](https://academic.oup.com/biomet/article/101/1/85/2365817) and graphical lasso for rank correlations [Nonparanormal](http://www.jmlr.org/papers/volume10/liu09a/liu09a.pdf). For simplicity, we only compare methods based on their ROC and AUC. Selection of tuning parameters is often criticized for being subjective.

 * Evaluation
 
      - Null model will be evaluated based on false positive rate of edge recovery (truth is no edge).
      
      - Alternative model will be evaluated by ROC and AUC.

