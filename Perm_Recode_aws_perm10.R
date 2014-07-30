print("R start time")
Sys.time()

library(reshape2)
library(dplyr)
library(plyr)
library(mail)

# multicore parallelization
library(doParallel)
#registerDoParallel(cores=2)

# read in phenotype data.
pheno <- read.table("2covarateN900V6_raw_Sean_Sex.txt", header=TRUE)

# Filter out subjects with same criteria as plink: remove LT_raw == -9 and ACE_4 == -9. Keep necessary columns
pheno <- subset(pheno, ACE_4 != -9 | LT_raw == -9, select = c('FID','LT_raw','ACE_4','SEX','Phenotype'))

# read in genetic data. Filtering out missing values based off Pheno filter above
gene <- read.table("Perm_Recode_Plink.raw", header=TRUE)
gene <- semi_join(gene, pheno, by='FID')

# Remove non-necessary columns, including or/not FID
drops.FID <- c('FID','IID','PAT','MAT','SEX','PHENOTYPE')
gene.noFID <- gene[, !names(gene) %in% drops.FID]
drops <- c('IID','PAT','MAT','SEX','PHENOTYPE')
gene.FID <- gene[, !names(gene) %in% drops]

# join genetics and phenotype together
# getting rid of unwanted columns for now (PCs, ChildAbuse, LT_trauma, etc)
#recode <- inner_join(pheno, gene) ==> seems to rearrange PAT, MAT columns. Revert back to plyr join
recode <- join(pheno, gene.FID, by="FID")

# melt the data for analysis
recode.melt <- melt(recode, id=c('FID','LT_raw','ACE_4','SEX','Phenotype'), 
                    variable.name="snps", value.name="snps_value")
#head(recode.melt)

# group by snps for bootstrapping
recode.melt.bysnps <- group_by(recode.melt, snps)

# analyzing the original hypothesis with main effects and interactions with lm()
# p values from this are a bit different than PLINK. Yet the orders are in the same range!!!
# !!!Check!!!
mod1 <- do(recode.melt.bysnps, failwith(NULL, lm), 
             formula= Phenotype ~ snps_value +SEX +LT_raw +ACE_4 +snps_value*LT_raw +snps_value*ACE_4)

# summarizing the model. correct one should be summary(), not anova(). see below in foreach() secition
mod1.summ <- lapply(mod1, summary)  # using summary gives you different results compared to anova()!! WHY??
mod1.summ.coef <- sapply(mod1.summ, coef)

# adding snps name to model1 
colnames(mod1.summ.coef) <- colnames(gene.noFID)

# extracting the p values for snpxLT_raw(27) and snpxACE4(28) as a new dataframe
origin.model.pval <- as.data.frame(mod1.summ.coef[27:28,])
origin.model.pval <- cbind(c("snpxLT_raw_pval", "snpxACE4_pval"), origin.model.pval)
colnames(origin.model.pval) <- c("gxe", colnames(mod1.summ.coef))

origin.model.pval.m <- melt(origin.model.pval, id='gxe', variable.name='snps', value.name='snps_pval')
origin.model.pval.long <- dcast(origin.model.pval.m, snps ~ gxe)

# > mod1.summ[1]
# [[1]]
# 
# Call:
#   f(formula = ..2, data = ..1)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -30.417  -7.378  -2.567   4.512  45.149 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       15.18680    1.64859   9.212  < 2e-16 ***
#   snps_value         1.52398    1.31135   1.162   0.2455    
# SEX                6.80922    1.18693   5.737  1.4e-08 ***
#   LT_raw             0.85689    0.07964  10.760  < 2e-16 ***
#   ACE_4              0.60440    0.23453   2.577   0.0102 *  
#   snps_value:LT_raw -0.01428    0.13685  -0.104   0.9169    
# snps_value:ACE_4   0.04344    0.37960   0.114   0.9089    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 11.56 on 746 degrees of freedom
# Multiple R-squared:  0.2053,  Adjusted R-squared:  0.1989 
# F-statistic: 32.11 on 6 and 746 DF,  p-value: < 2.2e-16

## === recruiting biglm() package to reduce memory size !!! Don't know how to extract coefficients though
# mod2 <- do(recode.melt.bysnps, failwith(NULL, biglm), 
#            formula= Phenotype ~ snps_value +SEX +LT_raw +ACE_4 +snps_value*LT_raw +snps_value*ACE_4)
# 
# mod2.list <- lapply(mod2, anova)    # anova() doesn't work with biglm
# mod2.summ <- lapply(mod2, summary)  # list view. summary works with biglm. needs more digits to display small p values

# > mod2.summ[1]
# [[1]]
# Large data regression model: f(...)
# Sample size =  753 
# Coef    (95%     CI)     SE      p
# (Intercept)       15.1868 11.8896 18.4840 1.6486 0.0000
# snps_value         1.5240 -1.0987  4.1467 1.3113 0.2452
# SEX                6.8092  4.4354  9.1831 1.1869 0.0000
# LT_raw             0.8569  0.6976  1.0162 0.0796 0.0000
# ACE_4              0.6044  0.1353  1.0735 0.2345 0.0100
# snps_value:LT_raw -0.0143 -0.2880  0.2594 0.1369 0.9169
# snps_value:ACE_4   0.0434 -0.7158  0.8026 0.3796 0.9089


## ============ using boot() bootstrap function with dplyr do() ==============


##### ========= trying out foreach() instead =======================
# Let's try the foreach() package to handle the do loop
#
library(foreach)

#library(iterators)
# can't iterate over factors!!!
# x <- foreach(recode.reduce.m$snps, .combine='c') %do% sum(recode.reduce.m$SEX)


##### ========= still using foreach() but with boot as iterator
# trying this from Stackoverflow
boot <- 5

boot_loop <- foreach(i=1:boot, .errorhandling='remove') %do% {
  # resampling phenotype columns
  # boot_pheno <- pheno[sample(nrow(pheno), replace=F),]
  
  # resampling only the Phenotype
  Phenotype <- pheno[sample(nrow(pheno), replace=F), 5]
  
  # blindly binding the resample phenotypes with original genotypes
  #boot_phenogeno <- cbind(boot_pheno, gene.noFID)       # don't like to use cbind(), still searching
  #boot_phenogeno_df <- tbl_df(boot_phenogeno)
  
  # resampling only the Phenotype
  boot_phenogeno <- cbind(pheno[,1:4], Phenotype, gene.noFID)
  
  # melt the genotype
  boot_phenogeno.m <- melt(boot_phenogeno, id=c('FID','LT_raw','ACE_4','SEX','Phenotype'), 
                           variable.name="snps", value.name="snps_value")
  # group_by the snps
  boot_phenogeno.m.bysnps <- group_by(boot_phenogeno.m, snps)
  
  # running the lm() thru do()
  boot_model <- do(boot_phenogeno.m.bysnps, failwith(NULL, lm), 
                   formula= Phenotype ~ snps_value +SEX +LT_raw +ACE_4 +snps_value*LT_raw +snps_value*ACE_4)
  
  # summarize bootstrap stats model
  boot_model.summ <- lapply(boot_model, summary.lm)   # same result as summary(). dumps into list
  boot_model.summ.coef <- sapply(boot_model.summ, coef) # sapply dumps into dataframe. Interested in row 27,28:: p val of gxe interactions
  # boot_model.summ[[1]]  # model verification purpose
  
  # finding the min p value across snps for snpxLT_raw(27) & snpxACE_4(28)
  boot.min.pval <- data.frame(min(boot_model.summ.coef[27,]), min(boot_model.summ.coef[28,]))
  
  # stacking these min p values up to test against original model
  #boot.min.pval.df <- NULL
  #boot.min.pval.df <- rbind(boot.min.pval.df, boot.min.pval)

  ## don't know how to query out p values with biglm()
#   boot_model.biglm <- do(boot_phenogeno.m.bysnps, failwith(NULL, biglm), 
#                          formula= Phenotype ~ snps_value +SEX +LT_raw +ACE_4 +snps_value*LT_raw +snps_value*ACE_4)
#   
#   boot_model.biglm.summ <- lapply(boot_model.biglm, summary)
  
  # Below seems to be wrong as this should be testing an aov() object, plus it's testing effects sequentially
  # Google: difference between summary.lm and summary.aov in r
  # http://stats.stackexchange.com/questions/28938/why-do-linear-regression-and-anova-give-different-p-value-in-case-of-consideri
  # http://stats.stackexchange.com/questions/20452/how-to-interpret-type-i-sequential-anova-and-manova/20455#20455
  # http://stats.stackexchange.com/questions/55571/r-differences-between-lm-and-aov
  # boot.model.summ.aov <- lapply(boot_model, summary.aov)  # same result as anova() 
  
  return(boot.min.pval)
}

print("foreach() end time")
Sys.time()

# stack list result into dataframe. ldply() works better than sapply(, cbind). stacks vertically.
boot.min.pval.df <- ldply(boot_loop)

# empirical p values calculation against original model
colnames(boot.min.pval.df) <- c("bs_snpxLT_raw_p", "bs_snpxACE4_p")

# using the ddply() way. parallel seems to be penalized due to small amount of data ~ 3500 rows
system.time(
  snpxLTraw <- ddply(filter(origin.model.pval.m, gxe=="snpxLT_raw_pval"), .(snps), mutate,
                     Emp_p_no = sum(abs(boot.min.pval.df$bs_snpxLT_raw_p) < abs(snps_pval)),
                     Emp_p = Emp_p_no / (boot + 1)
                     #, .parallel=TRUE
  )
)
                     
system.time(
  snpxACE4 <- ddply(filter(origin.model.pval.m, gxe=="snpxACE4_pval"), .(snps), mutate,
                    Emp_p_no = sum(abs(boot.min.pval.df$bs_snpxACE4_p) < abs(snps_pval)),
                    Emp_p = Emp_p_no / (boot + 1)
                    #, .parallel=TRUE
  )
)

write.csv(snpxLTraw, file="snpxLTraw_EMP_p.csv")
write.csv(snpxACE4, file="snpxACE4_EMP_p.csv")

# print out system time to see total calculation. also send out email for notification
print("R end time")
Sys.time()

sendmail("tehsheng@umich.edu", "R notice", "Calculation finished.\nFetch your data!")




# tried using dplyr chaining but didn't work
# origin.model.pval.m %.%
#   filter(gxe=="snpxLT_raw_pval") %.%
#   group_by(snps) %.%
#   summarise(
#     test = sum(abs(boot.min.pval.df$bootstrap_snpxLT_raw_p) < abs(snps_pval))
#     )

# tried using foreach() but decided to use ddply() instead
# test <- foreach()
# 
#   do(count, formula = abs(boot.min.pval.df$bootstrap_snpxLT_raw_p) < abs(snps_pval)) %.%
#   mutate(
#         Emp_p_snpxLT = 
#         )

# ====== using sample() directly without boot() help ==========
# will need to handle the parallelization later with foreach() or snow()
#



