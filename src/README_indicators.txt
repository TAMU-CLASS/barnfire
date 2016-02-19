This clustering module clusters indicators and uses labels to compute a generalized energy mesh:
    Create observations from subset of indicators, viz. indicators[strt:end]
    Cluster observations and extract labels
    Combine labels from each subset
    Merge labels into existing group structure
    Create generalized energy mesh using the dual of the input grid

NB. The default AgglomerativeClustering implementation is Ward. The example used on Scikit-learn's website is misleading because they have Ward and AgglomerativeClustering as two different examples, where AgglomerativeClustering uses an average linkage.
NB: Parallel K-means is broken on OS X. Always use option n_jobs=1 (default).
####
NB: Hierarchical agglomeration (Ward) without a specified connectivity has an O(N^3) complexity
and performs very slowly for large numbers of inputs. For v 0.15.1 ff (see below), the best way to
determine connectivity is to use kneighbors_graph with a small number of neighbors. In this use
case, clusters will only calculate distances to neighbor clusters, not all clusters. Tests were
run with numbers of neighbors up to 300 (1% of all points). Results indicate a linear scaling
of time-to-solution with the number of neighbors for both both kneighbors_graph and
AgglomerativeClustering. Tests with varying observation sizes show a O(n*N) complexity with n
the number of neighbors. This makes sense: there are N iterations to get to 1 cluster. The first
iteration requires looking at all of one's neighbors and computing distances. Subsequent iterations
only need to update the distances between the new cluster and its neighbors. The tests show that
the number of neighbors of a derived (new) cluster are approximately constant
(combining points doesn't increase the number of neighbors appreciably).
NB: If kneighbors_graph is used, the minimum number of neighbors required can be determined by
doing a connected components calculation. The neighbors graph should have one connected component.
NB: Changing the number of neighbors does have an effect on the final clustering, and may cause
one cluster to be split rather than another, leading to very different results. The agglomerative
clustering and kneighbors algorithms are deterministic (unlike the k-means implementation), so
keeping the same number of neighbors will produce consistent results. For most cases, the differences
in final clustering are small as a function of number of neighbors. However, consistency should be
maintained when varying the number of clusters.
NB: If rad_graph is used for connectivity, the radius used should be sufficiently large so the
number of connected components is 1. Tests indicate that the number of neighbors will not be constant
but will be around 25-30% on average of the total points. This does not avert O(N^3) behavior and
is slow in practice.
NB: Implementations of kneighbors_graph or Ward (now AgglomerativeClustering) in scikit-learn 0.14.1 were bad. The implementations in 0.15.1 are substantially improved. In 0.14.1, too few neighbors were chosen so cluster sizes of 1 occured when they shouldn't have.
NB: It seems as though a neighbor count between 20 and 30 works well, though a check for connected
###
components should be done (take max of 20 and minimum number of connected components, say).
NB. As of 0.16, Scikit-learn has a new Birch clusterer. Birch is similar to hierarchical agglomerative clustering with a connectivity matrix based on rad_graph, but Birch is built for speed and low memory usage. 
NB. One major disadvantage of this implementation of Birch is it requires the 'threshold' parameter to be set, which should be proportional to the observations' magnitude (not between 0 and 1). In fact, the threshold should be small enough so more clusters are produced by Birch than the final number of desired clusters. This means a nonlinear iteration on the threshold may be to be used. Since the relationship between threshold and final number of clusters is not smooth, bisection may need to be used.
NB. If the desire is to divide the observation input space such that each cluster has the same volume, the threshold should be proportional to the (V/N)^(1/d), where V is the volume filled densely with observations, N is the number of Birch clusters, and d is the dimensionality. This would split up the problem into equally-sizes pieces in observation space.
NB. Alternatively, MiniBatchKMeans or an connectivity-informed AgglomerativeClustering could be used to determine the initial clusters and their sizes (radii) could be computed. The smallest radius could be used as the input threshold to Birch.
NB. In practice, the final number of clusters produced by the Birch algorithm itself should be much larger than the number of desired clusters. Perhaps by a constant factor, say, 2-5. The last step of the sklearn implementation of Birch is to call AgglomerativeClustering to further cluster the output clusters from Birch into final clusters. This ought to be cheap because the final number of clusters should be reasonably small. The final number of clusters must be less than or equal to the number of output Birch clusters, which depends on the threshold (smaller threshold, more clusters).
NB. Birch does the following. First it builds a balanced tree to efficiently store the data. It keeps a triplet with the number of, linear sum of (a vector in features), and squared sum of (a norm over features) observations contained in the node for each node. Using these three pieces of information is sufficient to determine the similarity of a new data point to the existing data within a node (distance to centroids and variances). New data points are added one at a time. At each level of the tree, the triplet is used to determine to which child node a data point should go (the one with the smallest difference between centroids). The triplet for the node is updated as the data moves down the tree (the sum of two triplets is the triplet of the sums). Once the data point gets to a leaf, if the new data point is too far away from the data in the leaf, as determined by the threshold, a new leaf is created with the new data point. Otherwise, the data point is added to the data in the leaf. The square of the radius used for this comparison is the variance of the data in the leaf (computed in a clever way), normalized to the number of points in the leaf, but not normalized by the magnitude of the data itself (the threshold has the same units as the observations). As new leaves are added, if the number of children a node has increases beyond a branching_factor, the node is split into two nodes and the triplets for the new nodes are computed from the bottom up. Once all the data points are added, the tree is complete with some final number of leaves, each of which is a cluster of the data. These clusters are optionally condensed into a final number of clusters using another clustering algorithm.
NB. When nodes are rebalanced during insertion, one node will be split into two (potentially several times, going up the tree from the leaves). The two new nodes split the children. Which children go to which parent is determined in a not-straightforward way. The two most distant children are determined, based on centroids. Each other child goes with its closest child.
NB. Because of the way splitting is done, Birch depends on the input order of the data. This is true both with respect to the way the tree is (re)balanced and the way leaves are formed.
NB. Each node can only hold up to branching_ratio children, so each new data point is not compared to all possible nodes, but merely some of them. This can cause non-physical groupings to occur. To get around this problem, set branching_ratio to be large (larger than the final number of clusters).
NB. The Birch algorithm is fast because it avoids comparing the new data point to all existing data points/clusters (as occurs in the default AgglomerativeClustering implementation) or to all centroids (as occurs in the K-means implementation), but instead only compares to a logarithmically small number of nodes (which are like centroids). All expensive computations are done on a very local level, meaning data movement is minimized. For this reason, Birch can be applied to datasets that do not fit into RAM.
NB. It would be very straightforward to create a variant of Birch that used different criteria to decide which child of a node to take and how to divide children of a split node. For the former, it could be based on which child would get the smallest increase in total (not normalized to number of points) variance if the data point were to be added to it (currently, it's based on distance of centroids). For the latter, formal clustering (perhaps agglomerative on centroids) of the children could be done (currently, the two furthest centroids are used as team captains and children go to the parent corresponding to the closest captain's centroid, which does have the advantage of being L^infty-like). This could be used to spawn more than two parents, though the ideal number of parents to split into is likely problem-dependent.
NB. Birch can be used as a first step in determining the number of energy elements per coarse group. It can be used to estimate within-coarse-group variance more directly by estimating the required number of energy elements to get to a desired level of accuracy for that coarse group. To do this generally, determine the diameter of the data norm(max(obs,axis=0)-min(obs,axis=0)) and set the radius to be a small factor of this diameter. For this problem, knowing that observations are log10's of fluxes and cross sections, a threshold of 0.1 would create clusters whose effective radius would put two neighboring clusters a factor of 1.26 apart in flux/xs value, for example. Using this as the radius and n_clusters=None, the number of unique labels output by Birch would be proportional to the number of unknowns required to get to a desired accuracy level.
NB. Birch and hierarchical agglomerative clustering are quite comparable for a reasonably set threshold for Birch. For large enough threshold, Birch is much faster.








In old README file:

https://www-nds.iaea.org/wimsd/libraries.htm (WIMS energy group structures and other goodies)
NJOY: GROUPR: IGN option:
9 -- WIMS 69-group structure
10 - 187-group structure (good for fast or thermal)


The idea: you have a high-dimensional space (the UFG energy structure) that you want to represent with a lower dimensional space (the bands). You cluster the UFG into bands to minimize the variation in reproducing the angular flux at several points (spatial, angular) of interest. This is vector quantization (vq).
Use k-means clustering or agglomerative hierarchical clustering for vector quantization in Python.
There are many fast algorithms for k-means clustering, in particular.
For a C++ implementation of k-means, see http://www.mlpack.org/about.html.
Expectation maximization is another useful algorithm for clustering. See http://en.wikipedia.org/wiki/Cluster_analysis#Agglomerative_hierarchical_clustering.
PCA is (only) useful in decreasing the dimensional size (number of points). I have not found a good example of using PCA to determine clusters by itself. Instead, PCA can be used to (1) determine how many clusters are needed to explain some fraction of the variance or (2) to produce a smaller set on which to do k-means clustering.
