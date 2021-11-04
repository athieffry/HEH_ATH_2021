library(tidyverse)
library(magrittr) # for the pipes
library(reshape2) # for reshaping dataframes

'h' <- head # shortcut for the head() command

setwd('path/to/your/working/directory')



# 1. READ DATASET ####
# --------------------
TCs <- readxl::read_xlsx('CAGE_TC_counts.xlsx')

# a. look at wide format
TCs

# b. melt df into long format, using TC_id as "key", and plot
TCs %>%
    melt(id.vars='TC_id') %>%
    separate(col='variable', into=c('genotype', 'replicate'), sep='_', remove=F) %>%
    ggplot(aes(x=variable, y=value, fill=genotype)) +
           geom_boxplot() +
           scale_y_log10() +
           theme(axis.text.x=element_text(angle=45),
                 aspect.ratio=1) +
           labs(x='CAGE libraries', y='CAGE TC expression (counts)',
                title='This starts to look like something',
                subtitle='there is even place for a subtitle!',
                caption='Caption: Raw tag counts')



# 2. NORMALIZE EXPRESSION INTO TPM (Tags Per Million mapped tags) ####
# --------------------------------------------------------------------
# a. compute TPM factor for each library
TPM_factor <- TCs %>%
              select(-TC_id) %>%
              colSums() %>%
              divide_by(1000000)

# b. normalize expression to TPM using each library's TPM factor
TCs_TPM <- TCs %>%
           column_to_rownames('TC_id') %>%
           apply(1, function(x) divide_by(x, TPM_factor)) %>%
           t() %>%
           as.data.frame()

# c. Sanity check: all normalized libraries sum up to 1 million
TCs_TPM %>% colSums()

# d. Box plot of normalized CAGE expression
TCs_TPM %>%
    as.data.frame() %>%
    rownames_to_column('TC_id') %>%
    melt(id.vars='TC_id') %>%
    separate(col='variable', into=c('genotype', 'replicate'), sep='_', remove=F) %>%
    ggplot(aes(x=variable, y=value, fill=genotype)) +
           geom_boxplot() +
           scale_y_log10() +
           theme(axis.text.x=element_text(angle=45),
                 aspect.ratio=1) +
           labs(x='CAGE libraries', y='CAGE TC expression (TPM)',
                title='This starts to look like something',
                subtitle='there is even place for a subtitle!')



# 3. AVERAGE REPLICATES ####
# --------------------------
# a. how to group
TCs %>%
    melt(id_vars='TC_id') %>%
    separate('variable', c('genotype', 'rep'), sep='_') %>%
    group_by(TC_id, genotype)

# b. how to summarize across groups
TCs %>%
    melt(id_vars='TC_id') %>%
    separate('variable', c('genotype', 'rep'), sep='_') %>%
    group_by(TC_id, genotype) %>%
    summarize('mean_count'=mean(value))

# c. plotting mean expression across replicates using a density plot
    # for raw counts
    TCs %>%
        melt(id_vars='TC_id') %>%
        separate('variable', c('genotype', 'rep'), sep='_') %>%
        group_by(TC_id, genotype) %>%
        summarize('mean_count'=mean(value)) %>%
        ungroup() %>%
        ggplot(aes(x=mean_count, col=genotype)) +
               geom_density() +
               scale_x_log10() +
               labs(x='Counts (mean of replicates)',
                    title='Density plot')

    # for TPM-normalized counts
    TCs_TPM %>%
        rownames_to_column('TC_id') %>%
        melt(id_vars='TC_id') %>%
        separate('variable', c('genotype', 'rep'), sep='_') %>%
        group_by(TC_id, genotype) %>%
        summarize('mean_count'=mean(value)) %>%
        ungroup() %>%
        ggplot(aes(x=mean_count, col=genotype)) +
               geom_density() +
               scale_x_log10() +
               labs(x='TPM (mean of replicates)',
                    title='Density plot')

# d. how to ungroup and reset into wide format, in case it is is needed
TCs %>%
    melt(id_vars='TC_id') %>%
    separate('variable', c('genotype', 'rep'), sep='_') %>%
    group_by(TC_id, genotype) %>%
    summarize('mean_count'=mean(value)) %>%
    ungroup() %>%
    pivot_wider(names_from='genotype', values_from='mean_count')



# 4. MORE APPLY EXERCICES: PRINCIPAL COMPONENT ANALYSIS ####
# ----------------------------------------------------------
# a.
TCs_TPM %>% h

# b. remove observations (CAGE TCs) that have no expression variance across libraries
    # how many are we gonna discard?
    TCs_TPM %>% apply(1, function(x) var(x)!=0) %>% table()
    # get index of those to discard
    tokeep <- TCs_TPM %>% apply(1, function(x) var(x)!=0)
    # do the actual removal
    TCs_TPM_no0var <- TCs_TPM[tokeep, ]

# c. quick and dirty PCA in R
pca <- TCs_TPM_no0var %>%
       t() %>% # the expression matrix needs to be transposed before feeding to PCA function, see ?prcomp
       prcomp(scale.=T, center=T)

# d. look up the structure of the pca object
str(pca)

# e. PCA scatter plot
pca$x %>%
    as.data.frame() %>%
    rownames_to_column('sample') %>%
    separate('sample', c('genotype', 'rep'), sep='_') %>%
    ggplot(aes(x=PC1, y=PC2, col=genotype)) +
           geom_point(size=4)
