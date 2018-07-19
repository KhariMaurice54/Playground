# metabric_read.R

source('/Users/Roger/Box_Sync/DBMI/COSBBI-2018/Khari/K_play/PAM50/PAM50-R-code/pam50-centroids-read.R.R')
PAM50genes = rownames(PAM50centroids)
mb_gexprs = read.csv('/Users/Roger/Box_Sync/DBMI/COSBBI-2018/Khari/METABRIC/Metabric.gexprs.csv')
dim(mb_gexprs)
mb_gexprs[1:6,1:6]
mb_gexprs_names = names(mb_gexprs)
length(intersect(names(mb_gexprs), PAM50genes))  ## 45
setdiff(y=names(mb_gexprs), PAM50genes)  ## 5 missing
 ####    "CDCA1" "CXXC5" "KNTC2" "MIA"   "ORC6L"
p45 = intersect(names(mb_gexprs), PAM50genes)  ## 45
save(p45, file='Weakest_link/data/p45.rdata')

mb_gexprs_pts = mb_gexprs$X   # 1981
mb_gexprs_pam50 = mb_gexprs[
  c('X', 
    intersect(names(mb_gexprs), PAM50genes))] 

####  clinical data ####
dirMB = '/Weakest_link/misc/metabric 2/'
mb_clin = Metabric.clinical.features =
  read.csv(paste0(getwd(), dirMB,
       'Metabric.clinical.features.csv')
)
mb_clin_names = names(mb_clin)
mb_clin_pts = mb_clin$X    ## 1981
setdiff(y=mb_gexprs_pts, mb_clin_pts)   #### identical

#### 
mb_histo = Metabric.histo.features =
  read.csv(paste0(getwd(), dirMB,
                  'Metabric.PAM50.histo.csv'),
           header=FALSE
  )
names(mb_histo) = c('X', 'classif', 'histo')
mb_histo_names = names(mb_histo)
mb_histo_pts = mb_histo[[1]]    ## 1981
length(mb_histo_pts)
setdiff(y=mb_gexprs_pts, mb_histo_pts)   #### identical
table(mb_histo[c('classif', 'histo')])

mb_surv = Metabric.surv.features =
  read.csv(paste0(getwd(), dirMB,
                  'Metabric.survdata.DSS.csv'),
           header=TRUE
  )
names(mb_surv)[3] = 'cens'
require(survival)
kmfit = survfit(Surv(mb_surv$time, mb_surv$cens) ~ 1)
plot(kmfit, xlab='days', ylab='Survival')

mb = merge(mb_clin, mb_gexprs_pam50, by='X')
mb = merge(mb, mb_histo, by='X')
mb = merge(mb, mb_surv, by='X')
dim(mb)
save(mb, file='mb.rdata')
