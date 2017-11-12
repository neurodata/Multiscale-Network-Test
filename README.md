# Multiscale Network Test (MNT)

Testing independence between network toplogy and nodal attributes via diffusion maps applied to distance-based correlations.

## MNT Use

As an input of MNT, we require `igraph` object G, a vector or matrix representing nodal attributes X, types of test statistics, use of diffusion maps, a range of Markov time used for diffusion maps, and the number of permutation samples.

The following example tests independence between network topology of karate club and each member's faction. 

```
library(igraphdata)
data(karate)
G = karate
X = V(karate)$Faction
res <- NetworkTest.diffusion.stat(G, X, option = 1, diffusion = TRUE,
                                  t.range = c(0:5), n.perm = 500)
```
The above prints out a list containing all the mgc statistics under the observations and under 500 permuations at Markov time t=0,1,2,..,5. The next command chooses the optimal diffusion map embeddings with dafault of t=3 and then we print out the p-value of network dependence test and optimal t. 
```
res.optimal <- print.stat.optimal(res, default.t = 4)
print(res.optimal$pvalue)
print(res.optimal$alt.t)
```


## Reference

Details of methods used for network dependence testing in the `mgc` package are described in [arXiv](https://arxiv.org/abs/1703.10136).
