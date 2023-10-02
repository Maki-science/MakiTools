# MakiTools
Some useful tools I created for me and others that want to use them (at their own risk!).

I might advance this to a proper R package when I find some spare time. Some functions might also gain an own package once they are developed sufficiently.

So far the functions work for me as far as I tested. However, elaborated, or automated tests are not included yet.


# How to use it?
## Installation
You can install the (latest) development version of these tools from [GitHub](https://github.com/) with:

``` r
source("https://raw.githubusercontent.com/Maki-science/MakiTools/main/customFunctionsMaki.R")

```
All contained functions will be loaded into your environment, ready to use.
Simple tooltips are also available. However, additional tooltips, such as from help files are not included.
See the function annotations or file comments for further reading. I might expand these at some point if required.


## What it in there?
So far, there are several functions included. However, some already were outsourced to proper packages since they have been developed sufficiently and will soon be/have been published.

Here is a list of functions included:

  - eulerLotka() - calculates the population growth rate, including the first 3 reproduction cycles
  - DIndex.OneStep()[unfinished!] - thought to comprise all steps of the InARes framework (see my github repositories)
  - simD()[outdated!] - simulation function for performance simulation of the InARes framework (see my github repositories)
  - z_transform() - a group-wise z-transformation (standardization)
  - prep_data_dodged_bargraph() - restructures data for being used in a dodged bargraph (ggplot2). It was built when I required lots of it. However, there are also better ways (most times - but not always).
  - MakiCV() - performs a k-fold crossvalidation with n repetitions for (almost) all regression models. MRSE or Briar score are used according to model family (others or binomial). Additional parameters, such as number of iterations or seed, can be set.
  - get_thresholds()[under development!] - calculates thresholds (e.g., for food limitation), using binomal regression models
  - get_thresholds_gams()[under development!] - calculates threshold, using additive models


## Function use
Here I will provide more details on how to use the functions. I will expand the list as I have some spare time.





# Troubleshooting
If there are any issues or wishes for changes, you can open an issue here on Github (https://github.com/Maki-science/MakiTools/issues).


# Citation
To cite these functions in publications use:

Marvin Kiene.
MakiTools: an assembly of useful tools for statistics and ecology (2023).
https://raw.githubusercontent.com/Maki-science/MakiTools/main/customFunctionsMaki.R
