tree<-phytools:::lambdaTree(pbtree(n=10),lambda=0.5)
## create a translation table
## leaving a couple of single-taxon clades for fun
tip.label<-sample(tree$tip.label,8)
clade.label<-sample(c("Puerto Rico","Vieques","Virgin Islands"),8,
                    replace=TRUE)
N<-ceiling(runif(n=8,min=1,max=20))

N_sample <- runif(8, min= 2, max = 50)

test_tree <- star_tips(tree = woltree, clade.label = clade.label, nodes = nodes , tip.label = tip.label)
