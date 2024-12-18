# PhD projects code

Some code from the two main projects carried out thus far as part of my PhD

## Project 1 

Developed a metric using mean first passage time of a random walker to assess how invasible a landscape is. That is, how long will it take a random walker, starting at a source site, to traverse a network of patches and reach target patch for the first time. 

I try to predict the invasion time of a stochastic patch occupancy metapopulation model using General Additive Models.

I then compare this metric against a previously established metric known as Condatis. Condatis is based off an electirical circuit analogy.

This is developed to be used by conservationsist who need accurate and efficient metrics to "score" many thousands of alternative habitat network configurations. 

I used four algorithms to create a wide range of theoretical landscapes configurations so that our results would be robust against empirical networks. 

All code is in the programming language R version 4.1.2 (2021-11-01)

### Abstract:

Current protected habitat networks may be unsuitable to facilitate climate induced species range shifts. Designing permeable landscapes, through the creation or protection of habitat in the correct location, requires informed decision making. Computationally efficient metrics of landscape permeability enable conservationists to compare many thousands of alternative landscape configurations prior to assigning limited resources. Metapopulation models are considered biologically realistic, but their computational inefficiency makes them impractical to use for a large number of landscapes. “Condatis” has been developed as an efficient metric to bridge this gap. However, it is based on an analogy of a resistor network rather than any biological process. Here we propose an alternative metric, First Passage Time (FPT), which maintains the computational efficiency while portraying a more biologically meaningful process than Condatis. It is defined as the time taken for a random walker on a network of patches, starting at a source patch, to reach a target patch for the first time. Different statistical moments of FPT are then calculated. We compared the two metrics’ ability to predict (using statistical models) mean invasion time (as calculated by metapopulation models) on a set of theoretical landscape configurations. We found that while FPT statistics can be a good predictor of invasion time, Condatis was generally the strongest single predictor of invasion time. But combinations of FPT statistics could improve predictions significantly. Adding FPT statistics allowed predictions to be up to 10% more accurate than using Condatis alone. Since FPT statistics are inexpensive to compute, they could be a useful way to improve our assessments of landscape’s suitability to sustain range expansion under climate change. Our study highlights the application of First Passage Time as novel landscape summary metric, aiding practitioners to efficiently design landscapes robust to climate induced species range shifts. 

## Project 2 

Code which makes us chapter 2 and 3 of my thesis 

I build a simulation of a single population in discrete space (2-dimensional) and time (with a range of warming rates) with a latitudinal temperature gradient. 

Temporal and spatial environmental fluctuations are implemented with stationary Gaussian processes (with zero mean).

Dispersal is calculated by calculating the convolution of the population density lattice and a 2D negative exponential dispersal kernel using Fast Fourier Transformations

All code is in the python3 programming language 

### Abstract:

Temperature change can be highly volatile across space and time. Spatial variation creates climatic micro-refugia due to topographic complexity of landscapes, while temporal fluctuations (which are increasing with climate change) create opportunities for colonisation but also wipe out existing populations. A better understanding of this variation, and how spatial and temporal fluctuations interact, is essential to predict which species will be most impacted by current and future climate change. Using a single species population dynamics model in discrete time and space we find the emergence of a spatial lag between the latitudinal position of the populations range margin and the latitudinal point beyond which the temperature is too cold to sustain a reproducing population. We show how spatiotemporal fluctuations impact a species’ ability to track a shifting climate, through its effect on the size of the lag, along a latitudinal gradient at its colder range margin. 

