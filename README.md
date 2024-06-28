# MakiTools
Some useful tools I created for me and others that want to use them (at their own risk!).

Some of these functions may be published in a separate package at some point (when I have some time to build them).
Be aware, that the functions might be subject to further progress and development according to my needs or requests from other users. I'll do my best that all functionality will persit in future and just additional features will be added. However, in some cases, I might outsource functionality in a separate function if it seems appropriate at some point.

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


## What is in there?
So far, there are several functions included. However, some already were outsourced to proper packages since they have been developed sufficiently and will soon be/have been published.

Here is a list of functions included:

  - eulerLotka() - calculates the population growth rate, including the first 3 reproduction cycles (created for Daphnia population growth rates)
  - DIndex.OneStep()[unfinished!] - thought to comprise all steps of the InARes framework (see my github repositories)
  - simD()[outdated!] - simulation function for performance simulation of the InARes framework (see my github repositories)
  - z_transform() - a group-wise z-transformation (standardization)
  - prep_data_dodged_bargraph() - restructures data for being used in a dodged bargraph (ggplot2). It was built when I required lots of it. However, there are also better ways (most times - but not always).
  - MakiCV() - performs a k-fold crossvalidation with n repetitions for regression models based on base R and lme4 and gamm4 packages (lm, glm, lmer, glmer, gamm4). MRSE or Briar score are used according to model family (others or binomial). Additional parameters, such as number of iterations or seed, can be set.
  - MakiCV.nlme() - same as MakiCV(), but for models based on nlme and mgcv package (lme, gls, gam, gamm). It allows to include correlation or variance structures.
  - get_thresholds()[under development!] - calculates thresholds (e.g., for food limitation), using binomal regression models
  - get_thresholds_gams()[under development!] - calculates threshold, using additive models


## Function use
Here I will provide more details on how to use the functions. I will expand the list as I have some spare time.





# Troubleshooting
If there are any issues or wishes for changes, you can open an issue here on Github (https://github.com/Maki-science/MakiTools/issues).


# Citation
To cite these functions in publications use:

Marvin Kiene.
MakiTools: an assembly of useful tools for statistics and ecology (2024).
https://raw.githubusercontent.com/Maki-science/MakiTools/main/customFunctionsMaki.R
