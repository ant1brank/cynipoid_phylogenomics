library(ape)
library(cowplot)
library(tibble)
library(tidyr)
library(tidytree)
library(treeio)
library(ggplot2)
library(ggtree)

draw_tree <- function(tree, hexp, cyn_mono=TRUE, leg=TRUE, treescale_height=8.0)
{
    require(ggtree)
    require(ape)

    # Taxon names actually used in the trees
    taxonNames <- c(
    "Allox_ar","Phaen_vi","Callas_no","Gana_sp1","Lepto_bo","Lepto_cl",
    "Lepto_he","Parn_nig","Andr_cur","Andr_gro","Andr_qrm","Andr_qln",
    "Beloc_tr","Bior_pal","Calli_sp","Prot_spe","Cerop_ma","Irae_his",
    "Diast_ki","Peric_JH","Qwaq_sco","Syne_gif","Syne_jap","Syne_umb",
    "Syne_ito","Aula_tav","Isoc_cen","Aylax_hy","Hedic_le","Phana_JH",
    "Escha_ac","Dipl_spi","Pedia_ac","Cecin_ib","Nasonia","Orussus",
    "Micropl"
    )

    # Choose the corresponding display names we want to use
    displayNames <- c(
    "Alloxysta arc","Phaenoglyphis vil","Callaspidia not","Ganaspis sp","Leptopilina bou","Leptopilina cla",
    "Leptopilina het","Parnips nig","Andricus cur","Andricus gro","Andricus qrm","Andricus qln",
    "Belonocnema kin","Biorhiza pal","Neuroterus val","Protobalandricus spe","Ceroptres mas","Iraella his",
    "Diastrophus kin","Periclistus sp","Qwaqwaia sco","Synergus gif","Synergus jap","Synergus umb",
    "Synergus ito","Aulacidea tav","Isocolus cen","\"Aylax\" hyp","Hedickiana lev","Phanacis sp",
    "Eschatocerus aca","Diplolepis spi","Pediaspis ace","Cecinothofagus iba","Nasonia vit","Orussus abi",
    "Microplitis dem"
    )

    # Desired order of tips
    taxonOrder <- c(36,37,35,34,32,33,31,8,1,2,3,4,5,6,7,26,27,28,29,30,19,20,21,22,23,24,25,18,17,16,15,14,13,12,11,10,9)

    # Root the trees correctly
   tree <- ape::root(tree, outgroup="Orussus", edgelabel=TRUE)
   tree$node.label[1]<-""

    # Get indices for clades and branches we want to color
    tb <- as_tibble(tree)
    cynipidae <- MRCA(tb, "Andr_gro", "Escha_ac")$node
    figitidae <- MRCA(tb, "Parn_nig", "Lepto_he")$node
    diplolepidae <- MRCA(tb, "Dipl_spi", "Pedia_ac")$node
    paraulacidae <- tb$node[match("Cecin_ib",tb$label)]


    # Change display names and only show support < 100%
    for ( i in 1:length(tree$tip.label) )
     tree$tip.label[i] <- displayNames[ match(tree$tip.label[i],taxonNames) ]
     tree$node.label<-gsub("\\{","",tree$node.label)
    tree$node.label<-gsub("}","",tree$node.label)
    tree$node.label<-unlist(lapply(strsplit(tree$node.label,"/"),function(x){x[2]}))
    IC<-as.numeric(tree$node.label)
    tree$node.label=IC 
    # Fill empty values for tips
    IC[1]<-NA
    IC<-c(rep(NA,37),IC)
    if (leg==TRUE) {
        leg.pos <- c(.2,.78)
    }
    else {
        leg.pos <- "none"
    }

    font_size <- 3.0    # default is 3.88
    cyn_offset <- -0.14
    fig_offset <- -0.14
    dip_offset <- -0.18
    cyn_offset_text <- 0.02
    fig_offset_text <- 0.02
    dip_offset_text <- 0.02
    cyn_color <- 'darkseagreen2'
    # Draw tree
    ggtree(tree,ladderize = TRUE) +
        geom_tree(size=0.8) + geom_tiplab(aes(label = paste0("italic('", label, "')")), parse = TRUE, size = 2) +
	geom_treescale(x = 0.0, y = treescale_height, width = 0.1, fontsize=2.0) +
        theme(legend.position=c(x=.2, y=.7)) +
        scale_fill_gradient(low="#DC3220", high="#005AB5",na.value="NA")+
        geom_nodepoint(shape = 21, colour = "NA", aes(fill =IC) , size = 2.5, stroke = 0.5) +
        hexpand(hexp)
}
draw_tree2 <- function(tree, hexp, cyn_mono=TRUE, leg=TRUE, treescale_height=8.0)
{
    require(ggtree)
    require(ape)
    
    # Taxon names actually used in the trees
    taxonNames <- c(
    "Allox_ar","Phaen_vi","Callas_no","Gana_sp1","Lepto_bo","Lepto_cl",
    "Lepto_he","Parn_nig","Andr_cur","Andr_gro","Andr_qrm","Andr_qln",
    "Beloc_tr","Bior_pal","Calli_sp","Prot_spe","Cerop_ma","Irae_his",
    "Diast_ki","Peric_JH","Qwaq_sco","Syne_gif","Syne_jap","Syne_umb",
    "Syne_ito","Aula_tav","Isoc_cen","Aylax_hy","Hedic_le","Phana_JH",
    "Escha_ac","Dipl_spi","Pedia_ac","Cecin_ib","Nasonia","Orussus",
    "Micropl"
    )
    
    # Choose the corresponding display names we want to use
    displayNames <- c(
    "Alloxysta arc","Phaenoglyphis vil","Callaspidia not","Ganaspis sp","Leptopilina bou","Leptopilina cla",
    "Leptopilina het","Parnips nig","Andricus cur","Andricus gro","Andricus qrm","Andricus qln",
    "Belonocnema kin","Biorhiza pal","Neuroterus val","Protobalandricus spe","Ceroptres mas","Iraella his",
    "Diastrophus kin","Periclistus sp","Qwaqwaia sco","Synergus gif","Synergus jap","Synergus umb",
    "Synergus ito","Aulacidea tav","Isocolus cen","\"Aylax\" hyp","Hedickiana lev","Phanacis sp",
    "Eschatocerus aca","Diplolepis spi","Pediaspis ace","Cecinothofagus iba","Nasonia vit","Orussus abi",
    "Microplitis dem"
    )
    
    # Desired order of tips
    taxonOrder <- c(36,37,35,34,32,33,31,8,1,2,3,4,5,6,7,26,27,28,29,30,19,20,21,22,23,24,25,18,17,16,15,14,13,12,11,10,9)
    
    # Root the trees correctly
    tree <- ape::root(tree, outgroup="Orussus", edgelabel=TRUE)
    tree$node.label[1]<-""
    
    # Get indices for clades and branches we want to color
    tb <- as_tibble(tree)
    cynipidae <- MRCA(tb, "Andr_gro", "Escha_ac")$node
    figitidae <- MRCA(tb, "Parn_nig", "Lepto_he")$node
    diplolepidae <- MRCA(tb, "Dipl_spi", "Pedia_ac")$node
    paraulacidae <- tb$node[match("Cecin_ib",tb$label)]
    
    # Change display names and only show support < 100%
    for ( i in 1:length(tree$tip.label) )
    tree$tip.label[i] <- displayNames[ match(tree$tip.label[i],taxonNames) ]
    
    # get gene concordance factor
    gcf<-read.table("../iqtree/PerGeneAnalyses/GCF_iqtree/concordCynipoid.cf.stat",skip=17,header=T,fill=T)

    if (leg==TRUE) {
        leg.pos <- c(.2,.78)
    }
    else {
        leg.pos <- "none"
    }
    
    font_size <- 3.0    # default is 3.88
    cyn_offset <- -0.14
    fig_offset <- -0.14
    dip_offset <- -0.18
    cyn_offset_text <- 0.02
    fig_offset_text <- 0.02
    dip_offset_text <- 0.02
    cyn_color <- 'darkseagreen2'
    # fill empty values for tips
    gCF=c(rep(NA,38),gcf[,2])
    #draw tree
    ggtree(tree,ladderize = TRUE) +
    geom_tree(size=0.8) + geom_tiplab(aes(label = paste0("italic('", label, "')")), parse = TRUE, size = 2) +
    scale_fill_gradient(low="#DC3220", high="#005AB5",,na.value = "NA")+
    geom_nodepoint(shape = 21, aes(fill =gCF) , size = 2.5, stroke = 0.5,color="NA") +
    guides(color = guide_legend(override.aes = list(size = 4, shape = 15))) +
    theme(legend.position=c(x=.2, y=.7)) +
    hexpand(hexp)
}

t1 <- read.tree("../iqtree/PerGeneAnalyses/TC_IC_RAxML/RAxML_Corrected_Probabilistic_IC_Score_BranchLabelsCorrected.T5")

p1 <- draw_tree(t1, 0.22, TRUE, TRUE, 8.0)
t2<-read.tree("../iqtree/PerGeneAnalyses/GCF_iqtree/concordCynipoid.cf.branch")
p2 <- draw_tree2(t2, 0.22, TRUE, FALSE, 8.0)
labels <- c("A Internode Certainty", "B Gene Concordance Factor")
cowplot::plot_grid(p1,p2,ncol=2,labels=labels, label_size=10, hjust=0.0) + theme(plot.margin=unit(c(3,3,3,3), "pt"))


ggsave("Fig_S5.png", device="png")

