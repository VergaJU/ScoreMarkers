# ScoreMarkers

Function to score the expression of positive and negative markers in single cell RNA seq data.

The project aim is to score and label the cells from an anndata file based on the expression of positive and negative markers. The scoring function work as follow:<br/>
For each label the score is given by the formula below:<br/>
$$
\sum{(\alpha\times\frac{\sum_{i=1}^{n_p}{f(x_p)}}{n_p})+(\beta\times\frac{\sum_{i=1}^{n_n}{f(x_n)}}{n_n})}
$$

With:
- $\alpha$ = weight for positive markers
- $\beta$ = weight for negative markers
- $n_p$ = number of positive markers
- $n_n$ = number of negative markers
- $f{x_p}$ = function to score positive markers
- $f{x_n}$ = function to score negative markers


Where $f{x_p}$ is:

$$

$$