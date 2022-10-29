---
layout: post
title:  "Centrality Measures and The Networked SIR Model"
date:   2021-12-22 17:21:40 +0200
is_series: true
series_title: "SIR-SINDy"
---

<em>Note: In this, and all posts in this series, the terms "graph" and "network" are used interchangeably.
The same applies to the terms "node", "vertex" and "individual".</em>

In this post I will be talking about how observations of network centralities,
from simulations of the networked SIR model (previous post),
may be used as predictors when applying the SINDy systems identification method (another previous post).

Centrality measures have been considered before in the literature in connection with the networked SIR model,
for example in [Dekker, 2013],[Bucur et al., 2020] and [Yuan et al., 2013], but I have not seen papers where aggregates of centrality measures, for nodes in the different compartments, are tracked over the duration of a simulation,
which is what we will be looking at here.<br>

### What is a network centrality measure?
In the following let \\(V\\) be the vertex set of a finite graph \\(G\\), with cardinality \\(|V|=N\\).

Network centrality measures offer a way to quantify connectivity properties of a vertex \\(v_j \in V\\)
according to some criteria. Depending on the application, the measure establishes the importance of the
vertex from a connectivity standpoint.

Well known centrality measures include for example <em>degree centrality</em>, <em>betweenness centrality</em>,
<em>eigenvector centrality</em> and <em>closeness centrality</em> amongst others.

For example, the <em>degree centrality</em> \\(C_D(v_j)\\)of a vertex \\(v_j\\) is the number of edges that are adjacent to the it,
and so coincides with the vertex degree \\(deg(v_j)\\).

The betweenness centrality \\(C_B(v_j)\\) for \\(v_j \in V\\) is defined as the sum of the fractions of shortest paths (counting edges) in \\(G\\) going through
\\(v_j\\) between all vertice pairs \\(v_i,v_k \in V\\), where \\(i,k,j \in \{1,2,\ldots,N\}\\)
and \\(i\neq j,\, i < k \neq j\\), times a normalizing factor.

$$
C_B(v_j) = \frac{2}{(N-1)(N-2)}\, \sum_{i \neq j} \sum_{i < k \neq j} \frac{\sigma_{ik}(j)}{\sigma_{ik}}
$$

where \\(\sigma_{ik}(j)\\) denotes the number of shortest paths from \\(v_i\\) to \\(v_k\\) going through \\(v_j\\),
and \\(\sigma_{ik}\\) is the total number of shortest paths from \\(v_i\\) to \\(v_k\\) [Dekker, 2013].<br>
The normalizing factor \\(\frac{2}{(N-1)(N-2)}\\) above corresponds to the one
used in the <samp>igraph</samp> package, when computing the betweenness centrality in <samp>R</samp>
(using <samp>igraph::betweenness</samp> with the <samp>normalized=TRUE</samp> option),
but one may come across other factors in the literature.<br>
Betweenness centrality can for example be used to identify bottlenecks in traffic networks,
since traffic usually takes the shortest route from A to B.

### Aggregates of centralities
In my master's thesis I wanted to find some way to use dynamic centrality observations from simulations as predictors
when applying SINDy. My thinking was that an infected node's centrality measure, for example considering degree centrality,
would clearly influence the time derivative target that we want to approximate when using SINDy.
The solution I came up with, was to use the fraction of the centrality sum accounted for by the individuals
in the \\(S(t),\,I(t)\\) and \\(R(t)\\) vertex sets at time \\(t\\). If I may elaborate...

Let $$C_*(v_i)$$ be any per-node centrality measure for a graph \\(G\\) with \\(v_i \in V \\), then $$ C_*(G) = \sum_{j=1}^N C_*(v_j) $$ is a graph invariant. The quatity

$$
s_C(t) = \frac{1}{C_*(G)} \sum_{v_j \in S(t)} C_*(v_j)
$$

is an \\(s(t)\\)-like trajectory, but in terms of a fraction of a graph invariant, instead of fraction of a population.
The case being similar for the \\(i_C(t)\\) and \\(r_C(t)\\) fractions.
Depending on the measure \\(C_*\\) being used, the  different centrality-sum fractions (CSFs) above provide dynamical information about how the pathogen
in the networked SIR model is spreading through the population from a topological standpoint.
It is somewhat unclear exactly what information is being conveyed,
but the assumption is that the information is useful in a machine learning context, i.e using SINDy.

Another thing to consider, is that a model of the networked SIR dynamics using CSFs as predictors will have more than 3 equations &mdash; one additional for every CSF used as a predictor,
unless one assumes that the CSFs are known functions of time,
which would render the model less useful in an applied sense.<br>
The good news is, that the sum of CSFs is a conserved quantity, i.e. \\(s_C + i_C + r_C = 1\\),
which enables us to use som "tricks" when finding a model with SINDy.
More on this in the final post in this series.

### Simulation example
Let's see what these CSF trajectories look like when using degree centrality and betweenness centrality.<br>
I've made a small modification to the <samp>networked-sir.R</samp> code from 
<a href="{% post_url 2021-12-14-networked-sir %}">previous post</a>.
The <samp>networked_sir()</samp> function now takes an additional vector <samp>C</samp> as an argument which
gives some per-vertex centrality measure. The CSFs are then computed using the <samp>csf()</samp> function,
and returned in observation data frame.</p>

<samp>networked-sir-v2.R</samp>
```R
library(igraph)

csf <- function(X,C) {
  # Returns fractions of graph invariant sum(C)
  # in the S,I and R vertex sets
  # Parameters:
  # X : State vector of networked SIR model
  # C : Vertex centrality measures (same length as X)
  cg <- sum(C)
  sc <- sum(C[X == 0])/cg
  ic <- sum(C[X == 1])/cg
  rc <- 1-sc-ic
  return(list(sc = sc, ic = ic, rc = rc))
}

networked_sir <- function(X0, G, C, beta, gamma, tmax) {
  # Simulate the networked SIR model with arbitrary network G (igraph)
  # using the Gillespie stochastic simulation algorithm (SSA).
  # Parameters:
  # X0 : Initial configuration of state vector X0 in {0,1,2}^N,
  #   where {0,1,2} ~ {S,I,R} and N is the number of individuals
  # G : Network (igraph object) consisting of N nodes (individuals)
  # C : Vertex centrality measures
  # beta : per-edge infection rate
  # gamma : Recovery rate
  # tmax : Maximum time duration of simulation
  # Returns:
  # List with data frame [t,s(t),i(t),r(t),sc(t), ic(t), rc(t)], network G,
  # initial state vector X0, final state vector X, and parameters beta, gamma

  N <- length(X0)  # Number of individuals
  if(N != vcount(G) || N != length(C)) {
    error("Length of X0 has to equal length of C and number of vertices in G!")
  }

  lambda <- numeric(length = N)  # Allocate rate vector
  num_inf_neighbours <- numeric(length = N)  # Allocate number of infected neighbors vector

  # Set state X to initial state and time t to zero
  X <- X0;
  t <- 0

  # Keep track of number of individuals in each compartment.
  st <- length(which(X == 0))
  it <- length(which(X == 1))
  rt <- length(which(X == 2))

  # Compute CSFs
  csfs <- csf(X,C)

  # Initialize data matrix with observations and centrality fractions
  df <- data.frame(t = 0, s = st, i = it, r = rt,
                   sc = csfs$sc, ic = csfs$ic, rc = csfs$rc)

  dm <- as.matrix(df)  # Matrix is faster than data frame to update

  # Get adjacency list of G. This appears to be the fastest way
  # to find the neighbors of a particular individual
  adj_list <- igraph::as_adj_list(G)

  # Compute rates for all individuals
  for(i in 1:N) {
    if(X[i] == 0) {  # If susceptible
      neighbour_indices <- adj_list[[i]]
      num_inf_neighbours[i] <- length(which(X[neighbour_indices] == 1))
      lambda[i] <- beta*num_inf_neighbours[i]
    } else if(X[i] == 1) {  # If infected
      lambda[i] <- gamma
    } else {
      lambda[i] <- 0
    }
  }

  while(t < tmax) {
    lambda0 <- sum(lambda)
    if(lambda0 == 0) {  # If all rates are zero then we are done
      break
    }

    # Sample tau and event type
    tau <- rexp(1, rate = lambda0)
    t <- t + tau
    j <- sample(x = c(1:N), size = 1, prob = lambda/lambda0)

    neighbour_indices <- adj_list[[j]]  # Get neighbor indices of individual j

    # Infection event
    if(X[j] == 0) { # If individual j is susceptible
      X[j] <- 1  # Set individual j to infected
      lambda[j] <- gamma  # Update rate for individual j

      # Increment number of infected neighbors count for neighborhood
      num_inf_neighbours[neighbour_indices] <- num_inf_neighbours[neighbour_indices] + 1

      # Bookkeeping
      st <- st - 1
      it <- it + 1

    } else { # Recovery event
      X[j] <- 2  # Set individual j to recovered
      lambda[j] <- 0  # Set individual j rate to zero

      # Decrement number of infected neighbors count
      num_inf_neighbours[neighbour_indices] <- num_inf_neighbours[neighbour_indices] - 1

      # Bookkeeping
      it <- it - 1
      rt <- rt + 1
    }

    # Update rates for susceptible neighbours of individual j
    update_indices <- neighbour_indices[X[neighbour_indices] == 0]
    lambda[update_indices] <- beta*num_inf_neighbours[update_indices]

    # Update data matrix
    dm <- rbind(dm,c(t,st,it,rt,unlist(csf(X,C))))

  }  # end while

  # Return fractions instead of counts
  dm[,c("s","i","r")] <- dm[,c("s","i","r")]/N
  df <- as.data.frame(dm)
  row.names(df) <- NULL  # Remove row names
  return(list(df = df, G = G, C = C, X0 = X0, X = X, beta = beta, gamma = gamma ))
}
```

The following code simulates the same \\(s,i,r\\) trajectories (same random seed for the simulations),
but computes the CSFs using degree centrality and betweenness centrality respectively.

```R
  source("networked-sir-v2.R")

  # Population size
  N <- 1000

  repeat {
    # Degree sequence drawn from geometric distribution
    # where all vertices have at least one edge
    dseq <- rgeom(n= N, prob = 0.2) + 1
    if(is_graphical(dseq)) {
      break
    }
  }

  # Construct graph from degree sequence
  G <- degree.sequence.game(out.deg = dseq, method = "vl")
  # Centrality measures
  C1 <- degree(G)
  C2 <- betweenness(G)
  # Initial state vector
  X0 <- sample(c(0,1,2), N, replace = TRUE, prob=c(0.99,0.01,0))
  # Parameters
  beta <- 0.5
  gamma <- 1
  tmax <- 15

  set.seed(123)
  sim1 <- networked_sir(X0 = X0, G = G, C1, beta = beta, gamma = gamma, tmax = tmax)

  set.seed(123)
  sim2 <- networked_sir(X0 = X0, G = G, C2, beta = beta, gamma = gamma, tmax = tmax)

  dfm1 <- melt(sim1$df, id.vars = "t")
  dfm2 <- melt(sim2$df, id.vars = "t")

  plot1 <- ggplot(data = dfm1, mapping = aes(t, value, color=variable)) +
    geom_line() +
    theme_light()

  plot2 <- ggplot(data = dfm2, mapping = aes(t, value, color=variable)) +
    geom_line() +
    theme_light()
```

<figure>
  <img src="/static/images/net-sir-csf1.png" alt="simulation using degree centrality">
  <figcaption>Plot 1 of networked SIR model using degree centrality for CSFs.</figcaption>
</figure>

<figure>
  <img src="/static/images/net-sir-csf2.png" alt="simulation using degree centrality">
  <figcaption>Plot 2 of networked SIR model using betweenness centrality for CSFs.
  </figcaption>
</figure>

As can be seen from the plots above the fraction of the centrality-sum
accounted for by the nodes in the \\(S(t),I(t)\\) and \\(R(t)\\) vertex sets
are far from the same as the fraction of population in the same sets.<br>
Comparing \\(i(t)\\) and \\(i_C(t)\\) in the two plots, we see that around the peak of the epidemic
it is the nodes with high degree and betweenness that are infected.
Since an infected node with many edges (high degree) is likely to increase the infection rate of many susceptible neighbours,
and thereby likely the estimated derivative of \\(i(t)\\) and \\(s(t)\\), it seems natural to use the degree-CSFs as
predictors in regression where the target is the estimated derivative (i.e. the SINDy stage).<br>
The case is not so clear for betweenness-CSFs, since a vertex can have high betweenness but low degree.
Even so, CSFs using degree and betweenness centrality are strongly correlated ( \\( \approx 0.97 \\) for the two \\(i_C(t) \\) time-series above).

### Conclusion
It appears that using CSFs provides us with additional information about how a pathogen spreads through
a population in the networked SIR model.<br>
In the next post I will be using CSFs as predictors in sparse regression when using SINDy,
and leave it up to the LASSO algorithm to decide whether the CSF predictors should be included in the resulting model.<br>
It turns out that LASSO loves CSFs!


<h3>References</h3>
[Dekker, 2013] Dekker, A. H. Network centrality and super-spreaders in infectious disease epidemiology.
Proc. - 20th Int. Congr. Model. Simulation, MODSIM 2013 331–337 (2013) doi:10.36334/modsim.2013.a5.dekker.

[Bucur et al., 2020] Bucur, D. & Holme, P. Beyond ranking nodes: Predicting epidemic outbreak sizes by network centralities.
PLoS Comput. Biol. 16, 1–20 (2020).

[Yuan et al., 2013] Yuan, C., Chai, Y., Li, P. & Yang, Y. Impact of the network structure on transmission dynamics in complex networks. IFAC Proc. Vol. 13, 218–223 (2013).
