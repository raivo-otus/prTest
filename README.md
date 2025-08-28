
# prTest - Research compendium

<!-- badges: start -->
<!-- badges: end -->

This is a research compendium. This project investigates
the use of probabilisitc alternatives to classical tests
in the context of microbiome research in population study scale. 

The aims are to benchmark and test probabilistic alternative against classical
methods. Ensuring applicability to the scenarios that arise in the field. 

If found suitable, the end goal is to develop and package common workflows
into R/Bioconductor to help with widespread adoption and use. 

Raising awareness of the benefits and use of probabilistic methods can unlock 
powerful tools for researchers. 

# Problem space

Statistical hypothesis testing is a fundamental tool for assessing whether
observed differences between groups are likely to have arisen by chance.
In the univariate setting, classical approaches for two-sample and multi-group
comparisons—such as the t-test, Kruskal-Wallis, and ANOVA—are widely
used due to their simplicity, computational efficiency, and well-established
theoretical properties. However, these frequentist methods typically yield
point estimates and p-values, offering limited direct probabilistic
interpretation of parameter uncertainty.

In contrast, probabilistic (Bayesian) approaches treat unknown parameters as
random variables and yield full posterior distributions, enabling
richer inference such as direct probability statements about effect sizes,
estimation of uncertainty intervals, and incorporation of prior knowledge.
Despite these advantages, Bayesian methods can be more computationally
demanding and require careful modeling decisions.

This research focuses on comparing the performance, interpretability, and
practical trade-offs of probabilistic and classical methods for univariate
two-sample and multi-group hypothesis testing with microbiome data,
aiming to identify contexts in which one paradigm may offer clear benefits
over the other.

# Structure

Custom functions used in analysis reports reside in the R/ directory.

Analysis reports are found within the analysis/ directory.

Reports are rendered in github compatible markdown for ease of viewing.

# Common tests covered

Comparison of alpha diversity indices:

- Shannon diversity across two-groups and multiple-groups

