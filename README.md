Supporting Code and Data for "Sex-biased reduction in reproductive success drives selective constraint on human genes"

# Introduction

This repository contains the supporting code for our manuscript on the relationship between rare genetic variant burden as measured by s<sub>het</sub> and fertility. This repository consists of a few different resources necessary to replicate our findings:

1. An RStudio project consisting of three RMarkdown documents in `RMarkdown/`:
    + `SNVCalling_Filtering.Rmd` - Examples on how CNV QC and annotation was performed.
    + `CNVCalling_Filtering.Rmd` - Examples on how SNV annotation was performed.
    + `PhenotypeTesting.Rmd` - Code to replicate all main text figures, supplementary figures, and findings of the manuscript.
    
    These documents are easily loadable into RStudio by simply doing File -> Open Project and selecting the `UKBBFertility.Rproj` file. The first two documents (i.e. `*Calling_Filtering.Rmd` are **not intended** to be actually runnable, but are provided as examples of how we processed data as part of the project. `PhenotypeTesting.Rmd` on the other hand, is intended to be run, but requires the user to download UK Biobank participant protected data using their own [UK Biobank access](https://www.ukbiobank.ac.uk/register-apply/). While our manuscript is undergoing peer review, UK Biobank will not provide required files to run this document. We will update this document when data is made available thru data access with UK Biobank. We also provide in the directory `compiled_html/` html documents produced by [knitter](https://www.rforge.net/doc/packages/knitr/knit.html) which represent the data as we ran it on our system. Please view these documents in your browser if you want better formating with the following links:
    + [SNVCalling_Filtering](https://htmlpreview.github.io/?https://github.com/eugenegardner/UKBBFertility/blob/master/compiled_html/SNVCalling_Filtering.html)
    + [CNVCalling_Filtering](https://htmlpreview.github.io/?https://github.com/eugenegardner/UKBBFertility/blob/master/compiled_html/CNVCalling_Filtering.html)
    + [PhenotypeTesting](https://htmlpreview.github.io/?https://github.com/eugenegardner/UKBBFertility/blob/master/compiled_html/PhenotypeTesting.html)

2. Scripts used as part of Rmarkdown and other data processing in `scripts/`. Please see individual RMarkdown documents for more details.

3. Rawdata that is used as input for RMarkdown. This is provided as a tarball to save space. To use it please do the following:

    ```
    tar -zxf rawdata.tar.gz
    ```

    You should then be ready to use at least `PhenotypeTesting.Rmd` pending acquisition of UK Biobank data.

4. Java source code for tools that we created to do CNV and SNV/InDel annotation and QC in `src/`. This source code is provided with an [Eclipse IDE](https://www.eclipse.org/) project file to enable easy loading into eclipse. Both of these projects require external jars to compile:
    + Apache Commons [Math](http://commons.apache.org/proper/commons-math/), [CLI](http://commons.apache.org/proper/commons-cli/), and [Exec](http://commons.apache.org/proper/commons-exec/)
    + [htsjdk](https://github.com/samtools/htsjdk)

    Compiled jars which are runnable with a distribution of the **java14** JRE/JDK are also provided in `scripts/`. Please see CNV/SNV Calling and Filtering RMarkdowns for more information.

# Required Packages:

This project requires the following packages/dependencies:

* R:
    + biomaRt - Get gene lists we need (Need to install via [bioconductor](https://www.bioconductor.org/))
    + readxl - Read Supplemental Excel tables
    + data.table - Better than data.frame
    + patchwork - Arranging ggplots
    + broom - Makes getting odds ratios/betas/errors/p.values out of `glm()` much easier
    + meta - For doing meta analysis
    + mratios - Need this to calculate 95% CIs for ratios of two means
    + svglite - Need to create main text figures properly (ggsave doesn't like anything with an alpha shading)
    + tidyverse - Loads ggplot, tidyr, dplyr, and stringr
    + randomForest - R randomForest implementation
    + ROCR - Builds ROC curves
    + cvAUC - Determines CI for AUC from ROCR
    + gdata - Allows for accurate random sampling when making training sets for Random Forest
    + rcompanion - Used for calculating Naglekerke's pseudo r-squared for logistic regression
    + lubridate - For measuring time between patient birth/ICD-code incidence
* Perl:
    + None
* Java:
    + Apache Commons Math
    + Apache Commons CLI
    + Apache Commons Exec
    + [HTSJDK](https://github.com/samtools/htsjdk)
* Python
    + scipy
    + pandas

# Citing this Repo

If you use any of the materials in the repository, we would appreciate it if you cited our manuscript:

    (https://www.biorxiv.org/content/10.1101/2020.05.26.116111v2)
