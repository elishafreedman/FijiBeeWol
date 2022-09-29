tree<-phytools:::lambdaTree(pbtree(n=10),lambda=0.5)
## create a translation table
## leaving a couple of single-taxon clades for fun
tip.label<-sample(tree$tip.label,8)
clade.label<-sample(c("Puerto Rico","Vieques","Virgin Islands"),8,
                    replace=TRUE)
N<-ceiling(runif(n=8,min=1,max=20))

N_sample <- runif(8, min= 2, max = 50)

test_tree <- Backbone(tree = tree, clade.label = clade.label, nodes = N , tip.label = tip.label)


#manually setting depth
depths <- c(0.5*woltree$edge.length[170], 0.5*woltree$edge.length[149], 0.5*woltree$edge.length[143], 0.5*woltree$edge.length[135])
newtree <- phylo.toBackbone(woltree, trans = trans)

setEPS()
postscript("woltree.eps", width = 50, height = 100)
plot(woltree)
nodelabels(woltree$Nnode)
dev.off()

 star <- collapse.to.star(woltree, 164)
