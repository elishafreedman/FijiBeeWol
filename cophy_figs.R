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
make_phylo_class <- function(treedata = hosttree, sigfig = 2){

  d <- treedata@data[["posterior"]]
  tree <- list(edge =  treedata@phylo[["edge"]],
               Nnode = treedata@phylo[["Nnode"]],
               tip.label = treedata@phylo[["tip.label"]],
               edge.length = treedata@phylo[["edge.length"]],
               node.label = na.omit(round(
                 as.numeric(d), sigfig)))
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
  "/Users/freed/Documents/PhD_local_files/honours/cophylogeny/MCca_comb10%_fjWol_r.trees")
#woltree with specimen names
woltree_spec <- treeio::read.beast(
  "/Users/freed/Documents/PhD_local_files/honours/cophylogeny/MCca_comb10%_fjWol.trees")

woltree_spec@phylo[["tip.label"]]<- woltree_spec@phylo[["tip.label"]] %>%
  stringr::str_remove_all(., "wol") %>%
  stringr::str_remove_all(., "_") %>%
  stringr::str_replace_all(., "J0", "J") %>%
  stringr::str_replace_all(.,"fJ", "FJ") %>%
  stringr::str_remove_all(., "Consensus") %>%
  stringr::str_remove_all(., "consensus") %>%
  stringr::str_replace_all(., "cm", "CM_") %>%
  stringr::str_replace_all(., "H1051", "yCMR48_A06") %>%
  stringr::str_replace_all(., "HL054", "yCMR51_D06") %>%
  stringr::str_replace_all(., "HL455", "yCMR29_F03")

wol_phylo <- make_phylo_class(woltree_spec)

wol_nodes <- get_node_nums(wol_phylo, woltree_spec)

# hostree
hosttree <- treeio::read.beast(
  "/Users/freed/Documents/PhD_local_files/honours/cophylogeny/MCCA_EF1a2_COI2nd_nameChange.nex")

# rename tree
hosttree@phylo[["tip.label"]]<-
  hosttree@phylo[["tip.label"]] %>%
  stringr::str_replace_all(., "Homalictus_",
                           "Lasioglossum (Homalictus) ") %>%
  stringr::str_replace_all(., "sp_", "sp. ") %>%
  stringr::str_remove_all(., "_a")

host_phylo <- make_phylo_class(hosttree)
host_td<- drop.tip(object = hosttree,
                   tip = c("Lasioglossum (Homalictus) spNov_CdS_2",
                                        "Lasioglossum (Homalictus) ostridorsum_b",
                                        "Lasioglossum (Homalictus) tuiwawae_b" ,
                                        "Lasioglossum (Homalictus) tuiwawae_c",
                                        "Phantom_sp",
                                        "Lasioglossum (Homalictus) spNov_CdS_1",
                                        "Lasioglossum (Homalictus) fijiensis_b",
                                        "Homolictus_van_Levu_sp. CM25"),
                          trim.internal = TRUE,
                          rooted = FALSE,
                          collapse.singles = TRUE
                          )

host_phylo_td <- make_phylo_class(host_td)
host_nodes <- get_node_nums(tree = host_phylo_td, treedata = host_td)







# sampling data

coldata <- read.csv("/Users/freed/Documents/GitHub/Honours/Dorey_script/HomalictusCollectionData_2018.csv")

#### associations ####

assocs <- coldata[coldata$Specimen_code %in% wol_phylo$tip.label, c("Specimen_code", "Species_name")]

# clade conversion for species M

assocs$Species_name[which(assocs$Species_name
             == "Lasioglossum (Homalictus) sp. M")] <-
  "Lasioglossum (Homalictus) atritergus"
assocs$Species_name[which(assocs$Species_name
                          == "Lasioglossum (Homalictus) sp. R")] <-
  "Lasioglossum (Homalictus) concavus"
##### cophylogeny ####

cophy <- cophylo(wol_phylo, host_phylo_td,
                 assoc = assocs,
                 rotate = FALSE)
# quick check
setEPS()
pdf("cophy_rough.pdf", width = 50, height = 50)
plot(cophy, fsize = 3)
nodelabels.cophylo(which="left")
nodelabels.cophylo(which="right")
dev.off()
#  wol strains

# supergroup
cols <- colour_branches(tree = cophy,
                        col_l = c("orangered",
                                    "blue"),
                        col_r = NULL,
                        node_l_left = c(227, 129),
                        node_l_right = NULL)

# strains

links <- colour_links(tree = cophy,
                      names_from = cophy$assoc$Specimen_code,

                      names = c("CM_182",
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
pdf("cophy.pdf", width = 50, height = 80)
plot(cophy,
  lwd = 4,
  link.type = "curved",
  link.col = links,
  edge.col= cols,
  link.lty = "solid",
  link.lwd = 6,
  fsize = 3,
  use.edge.length = TRUE,
  align.tip.label = TRUE
)


# add posteriors

nodelabels.cophylo(cophy$trees[[1]]$node.label,
                   node = wol_nodes,
                     Ntip(cophy$trees[[1]]),
                   frame = "none",
                   adj = c(1, 1),
                   which = "left",
                   cex = 3
)


nodelabels.cophylo(cophy$trees[[2]]$node.label,
                    node = host_nodes,
                     Ntip(cophy$trees[[2]]),
                   frame = "none",
                   adj = c(1, 1),
                   which = "right",
                   cex = 3
)


dev.off()


##### smaller tree #####

# use same host tree

# tipdrop woltree

wol_drop <- wol_phylo$tip.label[!stringr::str_detect(wol_phylo$tip.label, "17FJ") &
                !stringr::str_detect(wol_phylo$tip.label, "CM_")]
wol_keep <-  c("17FJ127", "17FJ210", "17FJ98", "CM_182")


wol_drop_t <- wol_phylo$tip.label[!wol_phylo$tip.label %in% wol_keep &
                                    !wol_phylo$tip.label %in% wol_drop]


wol_td<- drop.tip(woltree_spec, wol_drop_t,
                   trim.internal = TRUE,
                   rooted = FALSE,
                   collapse.singles = TRUE
)

wol_phylo_td <- make_phylo_class(wol_td)
wol_nodes_col <- get_node_nums(tree = wol_phylo_td, treedata = wol_td)


# small tree associations

assoc_col <- assocs[assocs$Specimen_code %in% wol_keep,]

wHa_assocs <- data.frame(Specimen_code = rep("wHa", 10),
                         Species_name = c("Lasioglossum (Homalictus) hadrander",
                                          " Lasioglossum (Homalictus) kaicolo",
                                          "Lasioglossum (Homalictus) ostridorsum",
                                          "Lasioglossum (Homalictus) sp. S",
                                          "Lasioglossum (Homalictus) concavus",
                                          "Lasioglossum (Homalictus) groomi",
                                          "Lasioglossum (Homalictus) sp. F",
                                          "Lasioglossum (Homalictus) sp. J",
                                          "Lasioglossum (Homalictus) atritergus",
                                          "Lasioglossum (Homalictus) fijiensis"


                         ))
assoc_col <- rbind(assoc_col, wHa_assocs)

# change name of strains
wol_phylo_td$tip.label[wol_phylo_td$tip.label == "17FJ127"] <- "wFJA"
wol_phylo_td$tip.label[wol_phylo_td$tip.label == "17FJ210"] <- "wFJBa"
wol_phylo_td$tip.label[wol_phylo_td$tip.label == "17FJ98"] <- "wFJBb"
wol_phylo_td$tip.label[wol_phylo_td$tip.label == "CM_182"] <- "wFJBc"

# change name of strains
assoc_col[assoc_col$Specimen_code == "17FJ127",1] <- "wFJA"
assoc_col[assoc_col$Specimen_code == "17FJ210",1] <- "wFJBa"
assoc_col[assoc_col$Specimen_code == "17FJ98",1] <- "wFJBb"
assoc_col[assoc_col$Specimen_code == "CM_182",1] <- "wFJBc"

cophy_col <- cophylo(wol_phylo_td, host_phylo_td,
                 assoc = assoc_col,
                 rotate = FALSE)
# quick check
# quick check
setEPS()
pdf("cophy_col_rough.pdf", width = 50, height = 50)
plot(cophy_col, fsize = 3)
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
pdf("cophy_col.pdf", width = 40, height = 40)
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
