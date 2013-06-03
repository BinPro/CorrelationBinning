Timer unit: 1e-06 s

File: probin/binning/em.py
Function: cluster at line 8
Total time: 86.214 s

Line #      Hits         Time  Per Hit   % Time  Line Contents
==============================================================
     8                                           @profile
     9                                           def cluster(contigs, log_probability_func,fit_nonzero_parameters_func,cluster_count,centroids=None,max_iter=100, repeat=10,epsilon=1E-7):
    10         1           21     21.0      0.0      (max_clusters, max_clustering_prob,max_centroids) = (None, -np.inf, None)
    11        11          102      9.3      0.0      params = [(contigs, log_probability_func,fit_nonzero_parameters_func, cluster_count ,np.copy(centroids), max_iter,epsilon) for _ in xrange(repeat)]
    12                                               #pool = Pool(processes=cpu_count())
    13                                               #results = pool.map(_clustering_wrapper, params)
    14                                               #pool.close()
    15        11     86213906 7837627.8    100.0      results = [_clustering_wrapper(param) for param in params]
    16                                                   
    17         1           16     16.0      0.0      return max(results,key=lambda x: x[1])

File: probin/binning/em.py
Function: _clustering_wrapper at line 19
Total time: 86.2138 s

Line #      Hits         Time  Per Hit   % Time  Line Contents
==============================================================
    19                                           @profile
    20                                           def _clustering_wrapper(params):
    21        10     86213798 8621379.8    100.0      return _clustering(*params)

File: probin/binning/em.py
Function: _clustering at line 23
Total time: 86.2009 s

Line #      Hits         Time  Per Hit   % Time  Line Contents
==============================================================
    23                                           @profile
    24                                           def _clustering(contigs, log_probability_func,fit_nonzero_parameters_func, cluster_count ,p , max_iter, epsilon):
    25        10          223     22.3      0.0      if not np.any(p):    
    26        10     31566063 3156606.3     36.6          clustering,_, p = kmeans._clustering(contigs, log_probability_func, fit_nonzero_parameters_func, cluster_count ,p, max_iter=3,epsilon=epsilon)
    27        80          315      3.9      0.0          n = np.array([len(cluster) for cluster in clustering])
    28        10      4415330 441533.0      5.1          exp_log_qs, max_log_qs = _get_exp_log_qs(contigs,log_probability_func,p)
    29        10         1006    100.6      0.0          z = _expectation(contigs,log_probability_func,n,exp_log_qs)
    30        10      4400959 440095.9      5.1          prev_prob,_,_ = _evaluate_clustering(contigs, log_probability_func, p, z)
    31                                           
    32                                               else:
    33                                                   print >> sys.stderr, "Not implemented for EM to start with fixed p (centroids)"
    34                                                   sys.exit(-1)
    35                                                   
    36        10           22      2.2      0.0      prob_diff = np.inf
    37        10           21      2.1      0.0      prev_prob = -np.inf
    38        10           16      1.6      0.0      iteration = 0
    39                                               
    40        95          463      4.9      0.0      while(max_iter - iteration > 0 and prob_diff > epsilon):
    41        85         8367     98.4      0.0          z = _expectation(contigs,log_probability_func,n,exp_log_qs)
    42        85      8413608  98983.6      9.8          p = _maximization(contigs,fit_nonzero_parameters_func,z)        
    43        85     37374285 439697.5     43.4          curr_prob, exp_log_qs, max_log_qs = _evaluate_clustering(contigs,log_probability_func,p,z)        
    44        85         2735     32.2      0.0          n = np.sum(z,axis=0,keepdims=True)
    45                                                   
    46                                                   
    47        85          319      3.8      0.0          prob_diff = curr_prob - prev_prob
    48        85          141      1.7      0.0          (curr_prob,prev_prob) = (prev_prob,curr_prob)
    49        85          155      1.8      0.0          iteration += 1
    50                                               
    51                                               #Change back so curr_prob represents the highest probability
    52        10           18      1.8      0.0      (curr_prob,prev_prob) = (prev_prob,curr_prob)
    53        10          213     21.3      0.0      print >> sys.stderr, "iterations: {0}".format(iteration)
    54        10           43      4.3      0.0      if prob_diff < 0:
    55        10          206     20.6      0.0          print >> sys.stderr, "Previous probability {0} was greater than current probability {1}. Should never happen in EM".format(prev_prob,curr_prob)
    56        80          265      3.3      0.0      clustering = [set() for _ in xrange(len(p))]
    57        10          215     21.5      0.0      which_cluster = np.argmax(z,axis=1)
    58      3000         5871      2.0      0.0      for (contig,which) in izip(contigs,which_cluster):
    59      2990         9997      3.3      0.0          clustering[which].add(contig)
    60        10           14      1.4      0.0      return (clustering, curr_prob, p)

File: probin/binning/em.py
Function: _expectation at line 61
Total time: 0.008165 s

Line #      Hits         Time  Per Hit   % Time  Line Contents
==============================================================
    61                                           @profile
    62                                           def _expectation(contigs, log_probability_func, n, exp_log_qs):
    63                                               """
    64                                               Usage: _expectation(contigs,log_probability_func, p, n, exp_log_qs, max_log_qs)
    65                                               Return: The responsibility matrix for currenct p and n
    66                                               
    67                                               We are calculating the expression <z_{i,k}> = Q(theta_i|p_k)n_k / sum_{l=1}^K(Q(theta_i|p_l)*n_l)
    68                                                   where theta is the feature vector of contig_i, p_k is the feature vector of centroid_k and n_k is the expected number
    69                                                   of contigs in cluster k.
    70                                                   
    71                                                   the calculation of <z_{i,k}> can be expressed as:
    72                                                   
    73                                                   exp(log(Q(theta_i|p_k) - log(Q(theta_i|p_kmax))))*n_k / sum_{l=1}^K(exp(log(Q(theta_i|p_l))-log(Q(theta_i|p_max)))*n_l)
    74                                                   
    75                                                   We see that we get the exp_log_qs=exp(log(Q(theta_i|p_k))-log(Q(theta_i|p_max)) for all i and l and k 
    76                                                   and the max_log_qs = log(Q(theta_i|p_max) for all i for free from 
    77                                                   _evaluation_clustering from the previous iteration.
    78                                               """
    79                                               
    80        95         2133     22.5     26.1      exp_log_qs *= n
    81        95         6032     63.5     73.9      return exp_log_qs / np.sum(exp_log_qs,axis=1,keepdims=True)

File: probin/binning/em.py
Function: _maximization at line 83
Total time: 8.4126 s

Line #      Hits         Time  Per Hit   % Time  Line Contents
==============================================================
    83                                           @profile
    84                                           def _maximization(contigs, fit_nonzero_parameters_func, z):
    85                                               """
    86                                               Usage: _maximization(contigs, fit_nonzero_parameters_func, p, z)
    87                                               Return: The new centroids that maximize the probability for the current z values.
    88                                               
    89                                               We are calculating the expression p_{k,j} = sum_i(<z_{i,k}>*theta_{i,j}) / sum_j( sum_i(<z_{i,k}>*theta_{i,j}))
    90                                               
    91                                               """
    92        85      8412600  98971.8    100.0      return fit_nonzero_parameters_func(contigs,expected_clustering=z.T)

File: probin/binning/em.py
Function: _evaluate_clustering at line 95
Total time: 41.7739 s

Line #      Hits         Time  Per Hit   % Time  Line Contents
==============================================================
    95                                           @profile
    96                                           def _evaluate_clustering(contigs, log_probability_func, p, z):
    97                                               """
    98                                               Usage:  _evaluate_clustering(contigs,log_probability_func, p, z)
    99                                               Return: (clustering_prob, exp_log_qs, max_log_qs)
   100                                               
   101                                               calculate log L(theta|z,p) or log L(contigs| expected_clustering_freq,centroids)
   102                                               = sum_{i=1}^{N} (log (sum_{k=1}^{K} (<z_{i,k}>*Q(theta_{i}|p_{k}) ) ) )
   103                                               which translates into:
   104                                               sum_{i=1}^{N}( log( Q(theta_{i}|p_{kmax}) ) + log( sum_{k=1}^{K}( <z_{i,k}>*exp(log( Q(theta_{i}|p_{k}) ) - log(Q(theta_{i}|p_{kmax}) ) ) ) ) )    
   105                                               
   106                                               N is the number of contigs, K is the number of clusters.
   107                                               
   108                                               clustering_prob is the evaluation of the current clustering
   109                                               log_qs is an NxK matrix of all log (Q(theta_i|p_k)) values.
   110                                               max_log_qs is an Nx1 matrix where each row contains the max value over corresponding row in log_qs
   111                                           
   112                                               exp_log_qs and max_log_qs are used in the next iteration in _evaluation.     
   113                                               """
   114        95     41764083 439621.9    100.0      exp_log_qs, max_log_qs = _get_exp_log_qs(contigs,log_probability_func,p)
   115                                           
   116        95         9627    101.3      0.0      clustering_prob = np.sum(max_log_qs + np.log(np.sum(z*exp_log_qs)))
   117        95          147      1.5      0.0      return clustering_prob, exp_log_qs, max_log_qs

File: probin/binning/em.py
Function: _get_exp_log_qs at line 119
Total time: 46.0542 s

Line #      Hits         Time  Per Hit   % Time  Line Contents
==============================================================
   119                                           @profile
   120                                           def _get_exp_log_qs(contigs,log_probability_func,p):
   121       105          670      6.4      0.0      log_qs = np.zeros((len(contigs),len(p)))
   122     31500        86050      2.7      0.2      for i,contig in enumerate(contigs):
   123     31395     45938663   1463.2     99.7          log_qs[i] = np.fromiter((log_probability_func(contig,centroid) for centroid in p),dtype=float)
   124       105         3734     35.6      0.0      max_log_qs = np.max(log_qs,axis=1,keepdims=True)
   125       105        24911    237.2      0.1      exp_log_qs = np.exp(log_qs - max_log_qs)
   126       105          212      2.0      0.0      return exp_log_qs, max_log_qs

