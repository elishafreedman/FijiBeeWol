# This script was written by Elisha Freedman on the 20th of August 2023 to generate figures: 2 and S1
#for the manuscript titled "Widespread Wolbachia infection is correlated with increased mtDNA diversity in native bees across the Fijian archipelago."

# packages
library(phytools)
library(ape)
library(haplotypes)
library(pegas)
library(treeio)
library(tidytree)
library(phylotools)
library(TreeTools)
library(seqinr)
library(tidyverse)

##### functions ####
make_phylo_class <- function(treedata = hosttree, sigfig = 2, rev = FALSE){

  d <- treedata@data[["posterior"]]
  tree <- list(edge =  treedata@phylo[["edge"]],
               Nnode = treedata@phylo[["Nnode"]],
               tip.label = treedata@phylo[["tip.label"]],
               edge.length = treedata@phylo[["edge.length"]],
               edge.label = treedata@phylo[["edge.label"]],
               d = treedata@data[["posterior"]],
                 node.label = round(
                 as.numeric(d), sigfig))
  if(rev == TRUE){
    tree$node.label = rev(tree$node.label)
  } else{
    tree$node.label == tree$node.label
  }
  class(tree) <- "phylo"
  return(tree)
}


get_node_nums <- function(tree,
                          treedata){
  dataset <- treedata@data
  posts <- na.omit(data.frame(posterior = dataset$posterior,
                              node = dataset$node))
  nodes <- as.numeric(posts$node)
  return(nodes)
}
# colour tree with legends

colour_branches <-
  function(tree,
           col_l,
           col_r,
           node_l_left,
           node_l_right) {

    # tree 1

    side1 <- rep("black", nrow(tree$trees[[1]]$edge))
    if(is.null(node_l_left)){
      side1 <- side1
    }else{
    for (i in 1:length(node_l_left)) {
      nodes_l <- getDescendants(tree$trees[[1]], node_l_left[i])
      side1[sapply(nodes_l, function(x, y)
        which(y == x), y = tree$trees[[1]]$edge[, 2])] <-
        col_l[i]
    }
  }

    # tree 2

    side2 <- rep("black", nrow(tree$trees[[1]]$edge))
    if(is.null(node_l_right)){
      side2 <- side2
    } else {
    for(i in 1:length(node_l_right) == TRUE){
      nodes_r <- getDescendants(tree$trees[[2]], node_l_right[i])
      side2[sapply(nodes_r, function(x, y)
        which(y == x), y = tree$trees[[2]]$edge[, 2])] <-
        col_r[i]
    }
}
   col <- list(left = side1, right = side2)
   return(col)
  }


# colour links
colour_links <-
  function(tree,
           names_from,
           names,
           colours,
           main_col) {
    link_colours <- rep(main_col, nrow(tree$assoc))
    col_assign <-
      data.frame(cbind(link_colours, names_from))
    for (i in 1:length(names)) {
      col_assign[which(col_assign[,2] == names[i]), 1] <-
        colours[i]
    }
    col_assign <- c(col_assign[,1])
    return(col_assign)
  }

##### load files #####
# woltree without specimens
woltree_nspec <- treeio::read.beast(
   "/data/Fixed_nodePriors.tree")



woltree_nspec@phylo[["tip.label"]]<- woltree_nspec@phylo[["tip.label"]] %>%
  stringr::str_remove_all(., "wol") %>%
  stringr::str_remove_all(., "_") %>%
  stringr::str_remove_all(., " ") %>%
  stringr::str_replace_all(., "J0", "J") %>%
  stringr::str_replace_all(.,"fJ", "FJ") %>%
  stringr::str_remove_all(., "Consensus") %>%
  stringr::str_remove_all(., "consensus") %>%
  stringr::str_replace_all(., "cm", "CM_") %>%
  stringr::str_replace_all(., "H1051", "yCMR48_A06") %>%
  stringr::str_replace_all(., "HL054", "yCMR51_D06") %>%
  stringr::str_replace_all(., "HL455", "yCMR29_F03") %>%
  stringr::str_remove_all(., "reversed")

wol_phylo <- make_phylo_class(woltree_nspec)

wol_phylo[tip.label]


wol_nodes <- get_node_nums(wol_phylo, woltree_nspec)

# hostree
hosttree <- treeio::read.beast(
  "/data/CA10%_wbBees.tree")

hosttree@phylo[["tip.label"]]<-
  hosttree@phylo[["tip.label"]] %>%
  stringr::str_replace_all(., "Homalictus_",
                           "Lasioglossum (Homalictus) ") %>%
  stringr::str_replace_all(., "sp_", "sp. ") %>%
  stringr::str_remove_all(., "_a") %>%
  stringr::str_remove_all(., "reversed") %>%
  stringr::str_remove_all(., "_")

host_phylo <- make_phylo_class(hosttree)

##### smaller tree #####

#### small host tree
# tipdrop woltree

wol_drop <- wol_phylo$tip.label[!stringr::str_detect(wol_phylo$tip.label, "17FJ") &
                !stringr::str_detect(wol_phylo$tip.label, "CM_")]
wol_keep <-  c("17FJ127", "17FJ210", "17FJ98", "CM_182", "17FJ217", "17FJ213", "17FJ27")


wol_drop_t <- wol_phylo$tip.label[!wol_phylo$tip.label %in% wol_keep &
                                    !wol_phylo$tip.label %in% wol_drop]


wol_phylo_td<- drop.tip(wol_phylo, wol_drop_t,
                   trim.internal = TRUE,
                   rooted = FALSE,
                   collapse.singles = TRUE
)


# hosttree
host_drop <- host_phylo$tip.label[!stringr::str_detect(host_phylo$tip.label, "17FJ") &
                                   !stringr::str_detect(host_phylo$tip.label, "CM_")]
host_keep <-  c("17FJ127", "17FJ210", "17FJ98", "CM_182", "17FJ217", "17FJ213", "17FJ27")

host_drop_t <- host_phylo$tip.label[!host_phylo$tip.label %in% host_keep &
                                    !host_phylo$tip.label %in% host_drop]

host_phylo_td <- drop.tip(host_phylo, wol_drop_t,
                          trim.internal = TRUE,
                          rooted = FALSE,
                          collapse.singles = TRUE)
# small tree associations
coldata <- read.csv("/Users/freed/Documents/GitHub/Honours/Dorey_script/HomalictusCollectionData_2018.csv")
assocs <- coldata[coldata$Specimen_code %in% wol_phylo$tip.label, c("Specimen_code", "Species_name")]

cophy_small <- cophylo(wol_phylo_td, host_phylo_td,
                 assoc = assocs,
                 rotate = TRUE)
# quick check
setEPS()
pdf("/Users/freed/Documents/wol_paper versions/wol_paper_new_analysis/cophy_rough.pdf",
    width = 50, height = 50)
phytools::plot(cophy_small, fsize = 3)
nodelabels.cophylo(which="left")
nodelabels.cophylo(which="right")
dev.off()

### ready to plot

cols_col <- colour_branches(tree = cophy_col,
                        col_l = c("orangered",
                                  "blue"),
                        col_r = NULL,
                        node_l_left = c(65, 48),
                        node_l_right = NULL)

# strains

links_col <- colour_links(tree = cophy_col,
                      names_from = cophy_col$assoc$Specimen_code,

                      names = c("wFJBc",
                                "wFJBb",
                                "wFJBa",
                                "wFJA"
                      ),
                      colours = c("#882255",
                                  "forestgreen",
                                  "goldenrod1",
                                  "slateblue"),
                      main_col = "grey")

setEPS()
pdf("/Users/freed/Documents/wol_paper versions/wol_paper_new_analysis/cophy_col.pdf", width = 40, height = 40)
plot(cophy_col,
     lwd = 4,
     link.type = "curved",
     link.col = links_col,
     edge.col= cols_col,
     link.lty = "solid",
     link.lwd = 6,
     fsize = 3,
     use.edge.length = TRUE,
     align.tip.label = TRUE
)


# add posteriors

nodelabels.cophylo(cophy_col$trees[[1]]$node.label,
                   node = wol_nodes_col,
                   Ntip(cophy_col$trees[[1]]),
                   frame = "none",
                   adj = c(1, 1),
                   which = "left",
                   cex = 3
)


nodelabels.cophylo(cophy_col$trees[[2]]$node.label,
                   node = host_nodes,
                   Ntip(cophy_col$trees[[2]]),
                   frame = "none",
                   adj = c(1, 1),
                   which = "right",
                   cex = 3
)


dev.off()

#### large tree ###

host_phylo$tip.label <- host_phylo$tip.label %>%
  stringr::str_remove_all(., "reversed") %>%
  stringr::str_remove_all(., "_") %>%
  trimws()

woltree_nspec$tip.label <- woltree_nspec$tip.label %>%
  stringr::str_remove_all(., "reversed") %>%
  stringr::str_remove_all(., "_")

large_cophy<-cophylo(woltree_nspec, host_phylo,
                     rotate=TRUE)
setEPS()
pdf("/Users/freed/Documents/wol_paper versions/wol_paper_new_analysis/cophy_rough.pdf", width = 50, height = 50)
plot(large_cophy, fsize = 3)
nodelabels.cophylo(which="left")
nodelabels.cophylo(which="right")
dev.off()


cols <- colour_branches(tree = large_cophy,
                        col_l = c("orangered",
                                  "blue"),
                        col_r = NULL,
                        node_l_left = c(156, 129),
                        node_l_right = NULL)

links <- colour_links(tree = large_cophy,
                      names_from = large_cophy$assoc[,1],

                      names = c("CM182",
                                "17FJ98",
                                "17FJ210",
                                "17FJ127"
                      ),
                      colours = c("#882255",
                                  "forestgreen",
                                  "goldenrod1",
                                  "slateblue"),
                      main_col = "grey")

setEPS()
pdf("/Users/freed/Documents/wol_paper versions/wol_paper_new_analysis/large_cophy_curved.pdf", width = 50, height = 80)
plot(large_cophy,
     lwd = 4,
     link.type = "curved",
     link.col = links,
     edge.col= cols,
     link.lty = "solid",
     link.lwd = 6,
     fsize = 3,
     use.edge.length = TRUE,
     align.tip.label = TRUE,
     rotate = TRUE
)


nodelabels.cophylo(
  large_cophy$trees[[1]]$node.label[2:Nnode(large_cophy$trees[[1]])],
  2:Nnode(large_cophy$trees[[1]]) + Ntip(large_cophy$trees[[1]]),
  frame = "none",
  cex = 3,
  adj = c(1,1),
  which = "left"
)
nodelabels.cophylo(
  large_cophy$trees[[2]]$node.label[2:Nnode(large_cophy$trees[[2]])],
  2:Nnode(large_cophy$trees[[2]]) + Ntip(large_cophy$trees[[2]]),
  frame = "none",
  cex = 3,
  adj = c(1, 1),
  which = "right"
)

dev.off()

