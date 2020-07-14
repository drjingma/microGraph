# microGraph
Graphical models for microbiome data

 * Set the number of variables $p=200$ and vary the sample size $n=100,200,500$.
 
 * Run both null model (a. shuffle an available dataset, e.g. amgut1.filt from *SpiecEasi* package; b. generate from Dirichlet Multinomial) and parametric/alternative model (a. a log-normal distribution; b. generative copula model used in SPIEC-EASI). 
 
      - To simulate from Dirichlet Multinomial, we can first generate compositions from a log-normal distribution with mean $\mu$ and standard deviation 1.5: $\phi_{ij}\sim N(\mu_j, 1.5^2)$, where $\mu_j \sim {\rm Unif}[0,4]$. The compositions $$x_{ij} = \frac{\exp(\phi_{ij})}{\sum_{k=1}^p \exp(\phi_{ik})}.$$ Given the compositions, we can generate count data from a Dirichlet distribution. 
 
      - For log-normal distribution, you can use the same $\mu$ as above. For the inverse covariance matrix, generate an adjacency matrix from a random graph and entries following examples described in my code example. I can add the description of the details. 
      
      - For copula model, we can follow examples in SPIEC-EASI with the marginal being the zero-inflated negative binormial. Again this is also available in my code example.
 
 * Evaluation
 
      - Null model will be evaluated based on false positive rate of edge recovery (truth is no edge).
      
      - Alternative model will be evaluated by ROC and AUC.
      
      - To avoid being subjective, we do not plan to include comparisons for sensitivity and specificity at optimal tuning parameters. So no need to worry about stability selection. 

 * Implement methods including: Conet, SparCC, CCLasso, COAT, gCoda, SPIEC-EASI, SPRING. Note SPIEC-EASI is essentially running graphical lasso on CLR transformed data. Given a CLR transformed data, we can also apply other methods developed for Gaussian data, such as [SpaceJAM](https://academic.oup.com/biomet/article/101/1/85/2365817) and graphical lasso for rank correlations ([Nonparanormal](http://www.jmlr.org/papers/volume10/liu09a/liu09a.pdf)). 

 * Software tools
 
      - The package *Huge* is really useful. It has several built in functions for building ROC and calculating AUC. As long as we use the same sequence of tuning parameters, comparison across methods should be fair. 
      
      - I have included several code examples on running gCoda, Spiec-Easi, and SPRING.