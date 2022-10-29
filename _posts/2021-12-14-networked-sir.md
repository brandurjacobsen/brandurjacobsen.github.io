---
layout: post
title:  "The Networked SIR model"
date:   2021-12-14 17:21:40 +0200
is_series: true
series_title: "SIR-SINDy"
---

In the <a href="{% post_url 2021-12-09-sir-intro %}">previous post</a> we introduced the SIR model,
and argued that a networked model might be more appropriate if one seeks
to relax the assumption, that transmission of a pathogen can occur from any individual in the \\(I\\) compartment,
to any individual in the \\(S\\) compartment (with the compartments viewed as sets of individuals instead of counts),
and instead only is possible if there exists an edge, or symmetric relation between the individuals
in a network that models the social interactions of a population of size \\(N\\).

This of course introduces the new assumption, that the network is a useful model of the social interactions
(over the duration of the epidemic, in case of a static network, which we will assume the network is).

In this post I will show (with <samp>R</samp> code) how to perform stochastic simulations of a networked SIR model
using the Gillespie algorithm [Gillespie, 1977], while representing the network using the <samp>igraph</samp> package.

## Introduction
The networked SIR model is treated in e.g. [Kiss et al., 2017], and assumes that transmission
of a pathogen, from an infected to a connected susceptible individual,
takes place as a Poisson process with rate \\(\hat \beta\\) &mdash; the per-edge transmission rate,
while recovery of an infected individual takes place as a Poisson process with rate \\(\gamma\\),
independent of the network.

So, if a susceptible individual \\(s_j \in S, \; j \in \{1,2,\ldots,N\} \\) has \\(n_j \in \mathbb N,\; 0 \leq n_j \leq N-1\\)
infected neighbours, then transmission to \\(s_j\\) takes place as
a Poisson process with rate \\( \lambda_j = n_j\,\hat \beta\\).

A Poisson process is characterized by its exponentially distributed waiting times between events.<br>
Let \\( \tau_j \\) denote the waiting time until \\(s_j\\)
(with \\(n_j\\) infected neighbours) becomes infected, then \\(\tau_j \\) is exponentially distributed with density
$$
f_{\tau_j}(t;\,\lambda_j) = \lambda_j\,\exp^{-\lambda_j\,t} =  n_j\,\hat \beta \,\exp^{-n_j\,\hat \beta\,t}
$$

If at some time \\(t \geq 0\\) we have \\(s_j,s_k \in S \\), and \\(s_j,s_k\\) have infected neighbours.
Then which individual gets infected first?<br>
This is a case of so-called <em>racing exponentials</em>, and it turns out that the probability
\\( P(s_j\, \text{is infected first})=\frac{\lambda_j}{\lambda_j + \lambda_k} \\), i.e. proportional to the rate,
and similarly for \\(s_k\\). This fact is utilized in the Gillespie algorithm.

## The Gillespie algorithm
The Gillespie algorithm was initially invented in order to perform stochastic simulations of chemical reaction systems,
given by deterministic reaction-rate ODEs. Gillespie argued that stochastic effects could result in chemical reactions
having realizations that were far from what the deterministic models predicted, and thus advocated using stochastic simulations
as a complementary tool in conjunction with classical deterministic models when working with chemical reaction systems.<br>

We will not concern ourselves with these academic discussions, but rather focus on using the Gillespie algorithm to
simulate the networked SIR model.

We do this by simply considering infection and recovery events, instead of chemical reaction events,
since both types of systems assume exponentially distributed waiting times between events.
Then all that is required, is that we identify which types of events that can occur,
and their rate. 

Say we have \\(m\\) types of events. At the core of the Gillespie algorithm is the joint distribution
$$ P(\tau,j),\qquad \tau \in \mathbb R,\; j \in \{1,2,\ldots,m\} $$
which is the probability, given the state of the system at time \\(t\\),
that the next event will occur at time \\(t+\tau \\),
and will be an event of type \\(j\\).<br>
It turns out that we can sample from this distribution first by sampling \\(\tau\\)
from an exponential distribution with rate \\(\lambda\\) &mdash; the sum of all event rates at time \\(t\\),
i.e. \\( \lambda = \sum_{j=1}^m \lambda_j \\),
and then sample the type of event \\(j\\) with probability proportional the event rate \\(\lambda_j,\; j \in \{1,2,\ldots,m\}\\).<br>
But we are getting ahead of ourselves! We first need a network of individuals...

## Network representation
We can use the <samp>igraph</samp> package in <samp>R</samp> to represent the network of individuals in the networked SIR model.<br>
Every individual in the population of size \\(N\\) is represented by a vertex in the network.
The question is now: What is the topology of this network? Do we have to painstakingly construct this network by hand or what?<br>
Random graphs to the rescue! Unless we want to use a real-world network,
we simply sample the network uniformly from the set of connected networks with a
specified degree sequence, which as the term suggests, gives the degree (number of connections) of all the vertices.
We can do this with <samp>igraph</samp> package using the <samp>degree.sequence.game</samp> function. For example:

```R
library(igraph)

seq <- c(3,3,2,1,1) # Degree sequence
G <- degree.sequence.game(seq, method="vl")
plot(G)
```

<figure>
<img src="/static/images/network-example1.png" alt="Plot of example network" class="img-fluid">
<figcaption>Plot of an example network with prescribed degree sequence <samp>seq</samp>
produced by the code snippet above.</figcaption>
</figure>

Note that the degree sequence needs to be <em>graphical</em>, i.e. the sequence sum needs to be even,
and realizable as a simple graph &mdash; this can be determined using the
<a href="https://en.wikipedia.org/wiki/Havel%E2%80%93Hakimi_algorithm">Havel-Hakimi algorithm (Wikipedia)</a>.

So, if we have an assumption on the degree distribution of the population, from which the degree sequence can be sampled,
then we can, given that the degree sequence is realizable, get a <samp>igraph</samp> object \\(G\\) that we can use
in the implementation of the stochastic simulation.

## Simulation code
In the implementation below we have a state vector \\(X\\) of length \\(N\\).
The entries in the state vector can take on the values 0,1 or 2, corresponding
to a susceptible, infected or recovered individual.<br>
We also have a vector of rates \\(\lambda\\) of length \\(N\\), where
some entries may be zero, depending on the state of the corresponding individual,
and the state of the neighbours.<br>
After initialization we sample the joint distribution \\(P(\tau, j)\\) mentioned above,
and update \\(X\\) to reflect the new state of the system as a result of event \\(j\\).
We also update the rates that are affected by this new state,
i.e. some susceptible individuals may have more or fewer infected neighbours. Then repeat.<br>
The simulation stops and returns the observations
when there are no more infected individuals, or \\(t > tmax\\).

```R
  library(igraph)

  networked_sir <- function(X0, G, beta, gamma, tmax) {
    # Simulate the networked SIR model with arbitrary network G (igraph)
    # using the Gillespie stochastic simulation algorithm (SSA).
    # Parameters:
    # X0 : Initial configuration of state vector X0 in {0,1,2}^N,
    #   where {0,1,2} ~ {S,I,R} and N is the number of individuals
    # G : Network (igraph object) consisting of N nodes (individuals)
    # beta : per-edge infection rate
    # gamma : Recovery rate
    # tmax : Maximum time duration of simulation
    # Returns:
    # List with data frame [t,s(t),i(t),r(t)], network G,
    # initial state vector X0, final state vector X, and parameters beta, gamma

    N <- length(X0)  # Number of individuals
    if(N != vcount(G)) {
      error("Length of X0 has to equal number of vertices!")
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

    # Initialize data matrix with observations
    df <- data.frame(t = 0, s = st, i = it, r = rt)
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
      dm <- rbind(dm,c(t,st,it,rt))

    }  # end while

    # Return fractions instead of counts
    dm[,c("s","i","r")] <- dm[,c("s","i","r")]/N
    df <- as.data.frame(dm)
    row.names(df) <- NULL  # Remove row names
    return(list(df = df, G = G, X0 = X0, X = X, beta = beta, gamma = gamma ))
  }

```

Let's try out this code by sampling a random graph with 1000 vertices or individuals,
with a degree sequence drawn from a geometric distribution, with \\(\hat \beta = 0.5 \\)
and \\(\gamma = 1 \\). Initially 99% of the population is susceptible, while 1% is infected.

```R
library(igraph)
library(ggplot2)
library(reshape2)

source("networked-sir.R")

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

# Parameters
G <- degree.sequence.game(out.deg = dseq, method = "vl")
X0 <- sample(c(0,1,2), N, replace = TRUE, prob=c(0.99,0.01,0))
beta <- 0.5
gamma <- 1
tmax <- 15

sim <- networked_sir(X0 = X0, G = G, beta = beta, gamma = gamma, tmax = tmax)
dfm <- melt(sim$df, id.vars = "t")
ggplot(data = dfm, mapping = aes(t, value, color=variable)) +
  geom_line() +
  theme_light()
```

<figure>
<img src="/static/images/networked-sir-simulation.png" alt="Plot of an example simulation" class="img-fluid">
<figcaption>Plot of a simulation using <samp>networked-sir.R</samp> with prescribed degree sequence
drawn from a geometric distribution with \(p=0.2\) and theoretical mean \(\frac{1}{p}=5\).</figcaption>
</figure>

In the plot of the simulation, we see the stochastic effects in the \\(s,i\\) and \\(r\\) trajectories.
Compare this, for example, with the trajectories produced by the deterministc model in my
previous post: <a href="{% post_url 2021-12-09-sir-intro %}">The SIR Model</a>.

In a future post we will examine the effects that different degree distributions have on the dynamics of the networked SIR model,
and use the <a href="sindy-intro.html">SINDy framework</a> to obtain a parsimonious (deterministic and continuous)
differential equation model of the networked SIR dynamics given different degree distributions.

I hope you liked the post. Please leave any questions or comments in the commenting section below!

### References
[Kiss et al., 2017] Kiss, I. Z., Miller, J. C. & Simon, P. L. Mathematics of epidemics on networks:
From exact to approximate models. in Interdisciplinary Applied Mathematics vol. 46 1–413 (Springer, 2017).

[Gillespie, 1977] Gillespie, D. T. Exact stochastic simulation of coupled chemical reactions. J. Phys. Chem. 81, 2340–2361 (1977).
