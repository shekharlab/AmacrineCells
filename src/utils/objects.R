

ORTHOTYPES = c(
  "oAC1 [A2]", "oAC2 [VG3]", "oAC3", "oAC4", "oAC5", "oAC6 [A8]",
  "oAC7 [SEG]", "oAC8 [nGnG-gl]", "oAC9", "oAC10", "oAC11",
  "oAC12", "oAC13 [A17]", "oAC14 [A17]", "oAC15", "oAC16*", "oAC17*",
  "oAC18 [nNOS]", "oAC19 [CA2]", "oAC20 [NPY]", "oAC21*", "oAC22*",
  "oAC23*", "oAC24", "oAC25", "oAC26", "oAC27", "oAC28", "oAC29",
  "oAC30 [PDGFRA]", "oAC31", "oAC32", "oAC33 [VIP]", "oAC34", "oAC35",
  "oAC36* [nGnG-ga]", "oAC37* [CRH]", "oAC38* [nNOS]", "oAC39* [CA1]",
  "oAC40 [nNOS]", "oAC41", "oAC42* [SAC]"
)

dtColors = function(n = 435, omit = c('white', 'ivory', 'grey60', 'darkgrey', 'lightyellow', 'lightcyan', 'floralwhite', 'lightcyan1'), lighten.by = 0){
  colors = standardColors()[which(!standardColors() %in% omit)]
  if(length(colors) < n) n = length(colors)
  return(colorspace::lighten(colors[1:n], lighten.by))
}

labels2colors_dt = function(vector, ...){
  return(labels2colors(vector, colorSeq = dtColors(length(vector), ...)))
}

myPalette = function(n_colors, darken = 0.2){
  # Amacrine palette
  return(colorspace::darken(as.character(paletteer_c("grDevices::Spectral", n_colors), darken)))
}

AcPalette2 = function(n = 100, lighten.by = 0.4){
  colorspace::lighten(dtColors(n), lighten.by)
}

AcPalette3 = function(n, colors = c('deeppink', 'orange', 'chartreuse', 'cyan', '#584B9FFF', 'violet'), vary.by = 0.4){
  
  n.high = ceiling(n / length(colors))
  # n.low = floor(n / length(colors))
  col.list = unlist(lapply(colors, function(color) {
    colorRampPalette(c(color, colorspace::lighten(color, vary.by)))(n.high)
  }))
  
  col.list[1:n]
}

ClusterPalette2 = function(groups, colors = c('deeppink', 'orange', 'chartreuse', 'cyan', '#584B9FFF', 'violet', 'antiquewhite2'), 
                          debug = FALSE){
  
  levels = table(groups)[unique(groups)]
  col.list = unlist(lapply(seq_along(levels), function(i) {
    if(levels[[i]] > 1) {
      if(debug) browser()
      interweave(make_hcl_ramp(colors[[i]], levels[[i]]))
    } else {
      colors[[i]]
    }
  }))
  
  col.list
}


ClusterPalette = function(groups, colors = c('deeppink', 'orange', 'chartreuse', 'cyan', '#584B9FFF', 'violet', 'antiquewhite2'), 
                          vary.by = 0.2, dark.first = TRUE, debug = FALSE){
  
  levels = table(groups)[unique(groups)]
  col.list = unlist(lapply(seq_along(levels), function(i) {
    if(levels[[i]] > 1) {
      if(dark.first) {
        # if(levels[[i]] == 10 & debug) browser()
        # interweave(make_hcl_ramp(colors[[i]], levels[[i]]))
        interweave(colorRampPalette(c(colorspace::darken(colors[[i]], vary.by),
                                      # colors[[i]],
                                      colorspace::lighten(colors[[i]], vary.by)))(levels[[i]]))
      } else {
        interweave(colorRampPalette(c(colorspace::lighten(colors[[i]], vary.by), 
                                      # colors[[i]],
                                      colorspace::darken(colors[[i]], vary.by)))(levels[[i]]))
      }
    } else {
      colors[[i]]
    }
  }))
  
  col.list
}

DcPalette = function(n_colors, colors = c("violet", "deepskyblue", "chartreuse3", "orange",  "orangered","brown", 'grey')){
  # Double cone palette
  return(colorRampPalette(colors = colors)(n_colors))
}

# Global control
options(ggplot2.discrete.colour = lapply(c(4, seq(6,100, by = 4)), AcPalette3), 
        ggplot2.discrete.fill = lapply(c(4, seq(6,100, by = 4)), AcPalette3))

theme_set(
  theme_gray(base_family = "ArialMT")  # or any theme you like
)

species_palette = c(
  "Human" = 'deeppink2', 
  "Macaque" = 'deeppink4', 
  "Marmoset" = 'deeppink', 
  "MouseLemur"= 'deeppink3',  
  "TreeShrew"= 'darkred', 
  "Mouse"= 'darkorange',  
  "Rhabdomys"= 'darkorange3', 
  "Rat"= 'darkorange1',  
  "Peromyscus"= 'darkorange4', 
  "Squirrel"= 'goldenrod1', 
  "Cow"= 'cyan1',  
  "Sheep"= 'darkturquoise',  
  "Pig"= 'cyan2',  
  "Ferret"= 'cyan4', 
  "Opossum"= 'darkslategrey', 
  "Chicken"= 'chartreuse1', 
  "Lizard"= 'chartreuse3', 
  "Zebrafish"= 'antiquewhite3', 
  "Lamprey" = 'antiquewhite4'
)

species_palette2 = lighten(c(
  "Human" = 'deeppink3', 
  "Macaque" = 'deeppink1', 
  "Marmoset" = 'deeppink4', 
  "MouseLemur"= 'deeppink2',  
  "TreeShrew"= 'darkred', 
  "Mouse"= 'darkorange',  
  "Rhabdomys"= 'darkorange3', 
  "Rat"= 'darkorange1',  
  "Peromyscus"= 'darkorange4', 
  "Squirrel"= 'goldenrod1', 
  "Cow"= 'chartreuse2',  
  "Sheep"= 'chartreuse4',  
  "Pig"= 'chartreuse1',  
  "Ferret"= 'chartreuse3', 
  "Opossum"= 'darkslategrey', 
  "Chicken"= 'cyan1', 
  "Lizard"= 'cyan3', 
  "Axolotl" = 'blue',
  "Goldfish" = 'purple4',
  "Zebrafish"= 'darkviolet', 
  'Goldfish'= 'violet',
  "Killifish" = 'purple',
  'CatShark' = 'lightgrey',
  'Shark' = 'lightgrey',
  "Lamprey" = 'darkgrey'
), 0.15)


species_palette3 = lighten(c(
  "Human" = 'brown', 
  "Macaque" = lighten("deeppink1", amount = 0.3),
  "Marmoset" = 'deeppink4', 
  "MouseLemur"= 'deeppink2',  
  "TreeShrew"= 'darkred', 
  "Mouse"= 'darkorange',  
  "Rhabdomys"= 'darkorange4', 
  "Rat"= lighten('darkorange1', 0.3),
  "Peromyscus"= 'darkorange2', 
  "Squirrel"= 'goldenrod1', 
  "Cow"= 'chartreuse2',  
  "Sheep"= 'chartreuse4',  
  "Pig"= 'chartreuse1',  
  "Ferret"= 'chartreuse3', 
  "Opossum"= 'darkslategrey', 
  "Chicken"= 'cyan1', 
  "Lizard"= 'cyan4', 
  "Axolotl" = 'blue',
  "Newt" = 'dodgerblue',
  "Zebrafish"= 'darkviolet', 
  'Goldfish'= 'purple4',
  "Killifish" = 'violet',
  # 'CatShark' = 'lightgrey', # Don't use this one!
  'Shark' = 'lightgrey',
  "Lamprey" = 'darkgrey'
), 0.15)


setClass("SeuratAnnotation", 
         slots=list(cluster = "vector", 
                    count = "vector",
                    mean.feature = "vector",
                    mean.mt = "vector",
                    pct.doublet = "vector", 
                    classification = "vector",
                    lit.type = "vector", 
                    neuropeptides = "list",
                    proportion = "data.frame"
                    ))

SeuratAnnotation = function(SeuratObject, 
                            group.by = "seurat_clusters", 
                            peptide.list = as.list(NA)){
  
  metadata = Metadata(SeuratObject, group.by, "lit_type", "classification")
  
  cluster = metadata$seurat_clusters
  
  classification = metadata$classification
  
  pct.doublet = 1-(proportionTable(SeuratObject, 
                                feature.1 = "DF.classifications", 
                                feature.2 = group.by)[,"Singlet"])
  
  lit.type = metadata$lit_type
  
  # Extract cluster counts
  count = as.vector(table(SeuratObject[[group.by, drop = TRUE]]))
  
  # Extract mean nFeature
  mean.feature = MeanMetadata(SeuratObject, group.by, "nFeature_RNA")
  mean.mt = MeanMetadata(SeuratObject, group.by, "percent.mt")
  
  # Proportions is simply a cross-tabulation with each enrichment group
  enriched.table = as.data.frame.matrix(table(SeuratObject[[group.by, drop = TRUE]], SeuratObject[["enrichment", drop = TRUE]]))
  
  # Return object
  new("SeuratAnnotation", 
      cluster = cluster, 
      count = count,
      mean.feature = mean.feature,
      mean.mt = mean.mt,
      pct.doublet = pct.doublet, 
      classification = classification,
      lit.type = lit.type, 
      neuropeptides = peptide.list, 
      proportion = enriched.table/rowSums(enriched.table))
}

setMethod("show", "SeuratAnnotation", function(object) {
  
  # Collapse neuropeptides into a comma-delimited list
  neuropeptides = lapply(object@neuropeptides, function(x) {
      paste0(x, collapse = ", ") 
    }) %>% unlist
  
  df = data.frame(cluster = object@cluster, 
                  count = object@count,
                   mean.feature = object@mean.feature,
                   mean.mt = object@mean.mt,
                   pct.doublet = object@pct.doublet, 
                   classification = object@classification,
                   proportion = object@proportion,
                   lit.type = object@lit.type, 
                   neuropeptides = neuropeptides, 
                   check.names = FALSE)
  
  rownames(df) = NULL
  
  print(df)
})

# https://stackoverflow.com/questions/19612839/set-method-initialize-s4-class-vs-using-function
as.data.frame.SeuratAnnotation = function(object) {
  
  # Collapse neuropeptides into a comma-delimited list
  neuropeptides = lapply(object@neuropeptides, function(x) {
    paste0(x, collapse = ", ") 
  }) %>% unlist
  
  df = data.frame(cluster = object@cluster, 
                  count = object@count,
                  mean.feature = object@mean.feature,
                  mean.mt = object@mean.mt,
                  pct.doublet = object@pct.doublet, 
                  classification = object@classification,
                  proportion = object@proportion,
                  lit.type = object@lit.type, 
                  neuropeptides = neuropeptides, 
                  check.names = FALSE)
  
  rownames(df) = NULL
  
  return(df)
}

setClass("CelltypeAnnotation", 
         slots=list(gene = "vector", 
                    annotation = "vector", 
                    color = "vector"))

CelltypeAnnotation = function(gene, annotation, color, to.upper = FALSE, include.other = FALSE){
  
  # Save order if annotation is a factor
  if(is.factor(annotation)){
    order = levels(annotation)
    annotation = as.character(annotation)
  } else {
    order = unique(annotation)
  }
  
  # Add fake gene for other
  if(include.other){
    gene$Other = 'NOGENE'
    annotation = c(annotation, 'Other')
    order = c(order, 'Other')
    color = c(color, 'grey')
  }
  
  # Remove duplicate genes
  duplicates = getDuplicates(unlist(gene))
  gene = lapply(gene, function(x) x[!x %in% duplicates])
  
  # Check for duplicated genes 
  stopifnot(length(unique(unlist(gene))) == length(unlist(gene)))
  
  # Convert gene symbols to upper-case for seurat
  if(to.upper) gene = lapply(gene, toupper)
  
  # Replicate annotations by number of genes for each
  annotation = unlist(lapply(1:length(gene), function(i) {
    rep(annotation[i], length(gene[[i]]))
  }))
  
  # Replicate colors by number of genes for each
  color = unlist(lapply(1:length(gene), function(i) {
    rep(color[i], length(gene[[i]]))
  }))
  
  # Sort in order
  initial = data.frame(gene = unlist(gene), 
                       annotation = annotation, 
                       color = color)
  
  sorted = initial %>% arrange(factor(annotation, levels = order))
  
  # Return CelltypeAnnotation object
  new("CelltypeAnnotation", 
      gene = sorted$gene, 
      annotation = sorted$annotation,
      color = sorted$color)
}

setMethod("show", "CelltypeAnnotation", function(object) {
  print(data.frame(gene = object@gene, 
             annotation = object@annotation, 
             color = object@color))
})

as.data.frame.CelltypeAnnotation = function(object){
  data.frame(gene = object@gene, 
             annotation = object@annotation, 
             color = object@color)
}

subset.CelltypeAnnotation = function(object, genes){
  CelltypeAnnotation(gene = object@gene[match(genes, object@gene)], 
                     annotation = object@annotation[match(genes, object@gene)], 
                     color = object@color[match(genes, object@gene)])
}

# setGeneric("ColorCode", 
#            function(x, ..., verbose = TRUE) standardGeneric("myGeneric"),
#            signature = "x"
# )

ColorCode = function(object, vector){
  
  # Anything that is not annotated goes to NA
  vector[!vector %in% Annotation(object)] = NA
  
  for(i in 1:length(object@annotation)){
    vector[vector == object@annotation[i]] = object@color[i]
  }
  
  return(vector)
}

CelltypeOrder = function(object) unique(object@annotation)

setGeneric("Genes", function(x, ...) standardGeneric("Genes"))
setMethod("Genes", "CelltypeAnnotation", function(x, annotation = NULL) {
  if(!is.null(annotation)) x@gene[x@annotation == annotation] else x@gene
  })

setGeneric("Annotation", function(x) standardGeneric("Annotation"))
setMethod("Annotation", "CelltypeAnnotation", function(x) x@annotation)

setGeneric("Colors", function(x) standardGeneric("Colors"))
setMethod("Colors", "CelltypeAnnotation", function(x) x@color)

setGeneric("AnnotationColors", function(x) standardGeneric("AnnotationColors"))
setMethod("AnnotationColors", "CelltypeAnnotation", function(x) unique(x@color) %>% setNames(unique(x@annotation)))

# Amacrine subclass annotation
amacrine_annotation = CelltypeAnnotation(
                        gene = list(c("Chat", "Megf11", "Slc18a3", "Slc5a7", "Sox2"),
                                    c("NEFH", "Prkca", "Sdk1"), 
                                    c("NOS1"),
                                    c("Ddc", "Chl1", "Arhgdig"),
                                    c("TH"), 
                                    c("Slc17a7"), 
                                    c("Gjd2", "Prox1", "Dab1", "Nfia", "Dner", "Calb2"), 
                                    c("Satb2", "Ebf3"), 
                                    c("Slc17a8", "Sdk2")),
                          annotation = c("SAC",  
                                         "A17", 
                                         "nNOS", 
                                         "CA1",
                                         "CA2",
                                         "VG1", 
                                         "A2",  
                                         "SEG",  
                                         "VG3"), 
                        color = myPalette(9), 
                        to.upper = TRUE)
                          # color = c("red",
                          #                "orange", 
                          #                "gold", 
                          #                "chartreuse2", 
                          #                "darkgreen", 
                          #                "cyan", 
                          #                "skyblue", 
                          #                "darkblue", 
                          #                "magenta"))
                                         
amacrine_annotation2 = CelltypeAnnotation(
  gene = list(c("Chat", "Slc5a7", "Sox2"),
              c("NEFH", "Prkca", "Sdk1"), 
              c("NOS1"),
              c("Ddc", "TH"),
              c("TPBG"), 
              c("CARTPT"),
              c("SCGN"),
              c("SCG2"),
              c("PYGB"),
              c("CALB2", "npy", "cck", "calb1"),
              c("Tac1"),
              c("PDGFRA"),
              c("VIP"),
              c("Slc17a7"), 
              c("Gjd2", "Prox1", "Dab1", "Nfia", "Dner"), 
              c("SYT2"),
              c("Satb2", "Ebf3"), 
              c("Slc17a8", "Sdk2"), 
              c("PVALB"), 
              c("CRH"), 
              c("GBX2")),
  annotation = c("SAC",  
                 "A17", 
                 "nNOS", 
                 "CA1",
                 "CA2",
                 "CART",
                 "Spiny",
                 "Sec",
                 "Wiry",
                 "SemiL",
                 "SubP",
                 "T45",
                 "VIP",
                 "VG1", 
                 "A2",  
                 "A8",
                 "SEG",  
                 "VG3", 
                 "KT2", 
                 'CRH', 
                 'NNgaba'), 
  color = c("red",
                 "orange", 
                 "gold", 
                 "chartreuse2", 
                 "darkgreen", 
                 "grey30",
                 "grey80",
                 "grey30",
                 "grey80",
                 "grey30",
                 "grey80",
                 "grey30",
                 "darkred",
                 "cyan", 
                 "skyblue", 
                 "blue",
                 "darkblue", 
                 "magenta", 
                 "grey80", 
                  "grey80",
                  "grey30"))

amacrine_annotation_clean = CelltypeAnnotation(
  gene = list(c("CHAT", "SLC5A7"),
              c("PRKCA", "SDK1"), 
              c("NOS1"),
              c("DDC", "TH"),
              c("TPBG"), 
              c("PDGFRA"),
              c("VIP"),
              c("CRH"),
              c("GBX2"),
              c("NPY"),
              c("CARTPT", 'SLC35D3'),
              c('RXRG'),
              c('MAF'),
              
              c("GJD2"), 
              c("NFIB"),
              c("SATB2", 'EBF3'), 
              c("EBF2"), 
              c("SLC17A8", "SDK2"), 
              c("ROBO3"),
              c("TRHDE")),
  annotation = c("SAC",  
                 "A17", 
                 "nNOS", 
                 "CA1",
                 "CA2",
                 "PDGFRA",
                 "VIP",
                 "CRH",
                 "NNgaba",
                 "NPY",
                 "SLC35D3",
                 'RXRG',
                 'MAF',
                 
                 "A2",  
                 "A8",
                 "SEG", 
                 "NNgly",
                 "VG3", 
                 'ROBO3', 
                 'TRHDE'
                 ),
  color = c(head(myPalette(25), 13), tail(myPalette(25), 7))
)

# Major retinal class annotation
major_annotation = CelltypeAnnotation(
  gene = list(c("RBPMS", "SLC17A6", "NEFL", "NEFM"), #"POU6F2", "RBFOX3"), 
              c('VSX2', 'VSX1', 'OTX2', 'ISL1', 'GRM6', 'GRIK1'), #"PRKCA"), 
              c("PAX6", "TFAP2A", "TFAP2B", "TFAP2C", "GAD1", "GAD2", "SLC6A9"), 
              c("ONECUT1", "LHX1", "CALB1", "TPM3"), 
              c("PDE6H", "CRX", "ARR3", 'OPN1LW', 'OPN1MW', 'OPN1SW'), 
              c("SAG", "PDC", "RHO", 'NRL'), 
              # c("PDE6H", "CRX", "ARR3","SAG", "PDC", "RHO"),
              c("SLC1A3", "RLBP1",  "APOE", "APOB")), 
  annotation = factor(c("RGC",  # blue
                         "BC",  # chartreuse2
                         "AC", # cyan
                         "HC", # gold
                         "Cone",# orange
                         "Rod",  # red
                         # "PR", # red
                         "MG"), # magenta
                 levels = c("Rod", "Cone", "HC", "BC", "AC", "RGC", "MG")),  
  color = c("blue", 
                  "chartreuse2", 
                  "cyan", 
                  "gold", 
                  "red", 
                  "darkred", 
                  # "red",
                  "magenta"))

major_annotation_pr = CelltypeAnnotation(
  gene = list(c("RBPMS", 'RBPMS2', "SLC17A6", 'POU4F3', "NEFL", "NEFM"),
              c('VSX2', 'OTX2', 'ISL1', 'GRM6', 'GRIK1'), #"PRKCA"), 
              c("PAX6", "TFAP2A", "TFAP2B", "TFAP2C", "GAD1", "GAD2", "SLC6A9", 'SLC6A5'), 
              c("ONECUT1", "ONECUT2", "ONECUT3", "LHX1", "CALB1", "TPM3"), 
              c("PDE6H", "CRX", "ARR3","SAG", "PDC", "RHO"),
              c("SLC1A3", "RLBP1",  "APOE", "APOB")), 
  annotation = factor(c("RGC",  # blue
                        "BC",  # chartreuse2
                        "AC", # cyan
                        "HC", # gold
                        "PR", # red
                        "MG"), # magenta
                      levels = c("PR", "HC", "BC", "AC", "RGC", "MG")),  
  color = c("blue", 
            "chartreuse2", 
            "cyan", 
            "gold", 
            "red",
            "magenta"))
                  
major_annotation_custom = CelltypeAnnotation(
  gene = list(c("SAG", "PDC", "RHO"),
              c("PDE6H", "CRX", "ARR3"), 
              c("ONECUT1", "LHX1", "CALB1"), 
              c("VSX2", "CABP5", "GRIK1"), 
              c("TFAP2A", "GAD2", "SLC6A9"), 
              c("RBPMS", "SLC17A6", "NEFL"),  
              c("SLC1A3", "RLBP1",  "APOE"), 
              c("Other")), 
  annotation = factor(c("Rod", "Cone", "HC", "BC", "AC", "RGC", "MG", "Other"),
                      levels = c("Rod", "Cone", "HC", "BC", "AC", "RGC", "MG", "Other")),  
  color = c(myPalette(7), "grey"))

Rod_markers = Genes(major_annotation, annotation = "Rod")
Cone_markers = Genes(major_annotation, annotation = "Cone")
HC_markers = Genes(major_annotation, annotation = "HC")
BC_markers = Genes(major_annotation, annotation = "BC")
AC_markers = Genes(major_annotation, annotation = "AC")
RGC_markers = Genes(major_annotation, annotation = "RGC")
MG_markers = Genes(major_annotation, annotation = "MG")

phylogenetic_order = c("Human","Macaque","Marmoset","MouseLemur","TreeShrew",
                       "Mouse","Rhabdomys", "Rat","Peromyscus","Squirrel",
                       "Ferret","Pig","Cow","Sheep","Opossum","Chicken",
                       "Lizard", "Zebrafish", "Goldfish", "Lamprey")

phylogenetic_order2 = c("Human","Macaque","Marmoset","MouseLemur","TreeShrew",
                       "Mouse","Rhabdomys", "Rat","Peromyscus","Squirrel",
                       "Cow","Sheep","Pig","Ferret","Opossum","Chicken",
                       "Lizard", 'Axolotl', 'Newt', "Zebrafish", "Goldfish", 
                       'Killifish', 'Shark', 'CatShark', "Lamprey")

phylogenetic_order_pretty = c("Human","Macaque","Marmoset","Mouse lemur","Tree shrew",
                        "Mouse","Rhabdomys", "Rat","Peromyscus","Squirrel",
                        "Cow","Sheep","Pig","Ferret","Opossum","Chicken",
                        "Lizard", 'Axolotl', 'Newt', "Zebrafish", "Goldfish", 
                        'Killifish', 'Catshark', "Lamprey")

SPECIESFILES = c(Human = "Human", 
                 Macaque = "Macaque", 
                 Marmoset = "Marmoset", 
                 MouseLemur = "MouseLemur", 
                 TreeShrew = "TreeShrew", 
                 Mouse = "Mouse",
                 Rhabdomys = "Rhabdomys", 
                 Rat = "Rat", 
                 Peromyscus = "Peromyscus", 
                 Squirrel = "Squirrel", 
                 Cow = "Cow", 
                 Sheep = "Sheep", 
                 Pig = "Pig", 
                 Ferret = "Ferret", 
                 Opossum = "Opossum", 
                 Chicken = "Chicken_reclustered", 
                 Lizard = "Lizard_ncbi", 
                 Axolotl = 'Axolotl',
                 Newt = 'Newt',
                 Zebrafish = "Zebrafish", 
                 Goldfish = "Goldfish", 
                 Killifish = "Killifish_ncbi", 
                 CatShark = "CatShark", 
                 Lamprey = "Lamprey")

species.names = names(SPECIESFILES)
species.names.pretty = species_names_pretty(species.names)

FINALFILES = lapply(c(Human = "HumanAC_v6.rds", 
                      Macaque = "MacaqueAC_v6.rds", 
                      Marmoset = "MarmosetAC_v6.rds", 
                      MouseLemur = "MouseLemurAC_v6.rds", 
                      TreeShrew = "TreeShrewAC_v6.rds", 
                      Mouse = "MouseAC_v6.rds",
                      Rhabdomys = "RhabdomysAC_v6.rds", 
                      Rat = "RatAC_v6.rds", 
                      Peromyscus = "PeromyscusAC_v6.rds", 
                      Squirrel = "SquirrelAC_v6.rds", 
                      Cow = "CowAC_v6.rds", 
                      Sheep = "SheepAC_v6.rds", 
                      Pig = "PigAC_v6.rds", 
                      Ferret = "FerretAC_v6.rds", 
                      Opossum = "OpossumAC_v6.rds", 
                      Chicken = "ChickenAC_v6.rds", 
                      Lizard = "LizardAC_no_HCs_v6.rds", 
                      Axolotl = "AxolotlAC_v6.rds",
                      Newt = "NewtAC_v6.rds",
                      Zebrafish = "ZebrafishAC_no_HCs_v6_converted.rds", 
                      Goldfish = "GoldfishAC_v6.rds", 
                      Killifish = "Killifish_ncbiAC_v6.rds",
                      CatShark = "SharkAC_v6.rds",
                      Lamprey = "LampreyAC_v6.rds"), 
                    function(file) paste0('../../Species_Objects/', file))

# Photoreceptor type annotation
pr_annotation = CelltypeAnnotation(
  gene = list(c("full"), 
              c("opn1mw4/opn1lw1"),
              c("ZEB2", "zeb2a", "zeb2b"),
              c("opn1sw1", "OPN1SW"), 
              c("opn1sw2", "OPN2SW", "LOC132767847"),
              c("opn1mw1", "opn1mw2", "opn1mw3", "opn1mw4", "OPN1MSW", "RH2", "OPN1MW", "LOC132773706"),
              c("opn1lw1", "opn1lw2", "OPN1LW", "LOC132767849", "ENSSTOG00000024701"),
              c("THRB", "thrb"),
              # c("principal"),
              c("MYLK", "mylk", "STBD1", "stbd1"),
              c("RHO", "rhol")),
  annotation = factor(c("full",
                        "opn1mw4/\nopn1lw1",
                        "singleCone",
                        "UV", 
                        "blue", 
                        "green", 
                        "red", 
                        "principal",
                        # 'principle',
                        "accessory",
                        'rod'), 
                      levels = rev(c("rod", "singleCone", "UV", "blue", "green", "red", "opn1mw4/\nopn1lw1", "principal", "accessory", "full"))),  
  color = c("tan", "gold", "white", "violet", "deepskyblue", "chartreuse3", "orangered", "brown", "orange", "lightgrey"), 
  to.upper = FALSE)

# With calbindin
pr_annotation2 = CelltypeAnnotation(
  gene = list(c("full"), 
              c("opn1mw4/opn1lw1"),
              c("ZEB2", "zeb2a", "zeb2b"),
              c("opn1sw1", "OPN1SW"), 
              c("opn1sw2", "OPN2SW", "LOC132767847"),
              c("opn1mw1", "opn1mw2", "opn1mw3", "opn1mw4", "OPN1MSW", "RH2", "OPN1MW", "LOC132773706"),
              c("opn1lw1", "opn1lw2", "OPN1LW", "LOC132767849", "ENSSTOG00000024701"),
              c('CALB1', 'calb1', "THRB", "thrb"),
              # c("principal"),
              c("MYLK", "mylk", "STBD1", "stbd1"),
              c("RHO", "rhol")),
  annotation = factor(c("full",
                        "opn1mw4/\nopn1lw1",
                        "singleCone",
                        "UV", 
                        "blue", 
                        "green", 
                        "red", 
                        "principal",
                        # 'principle',
                        "accessory",
                        'rod'), 
                      levels = rev(c("rod", "singleCone", "UV", "blue", "green", "red", "opn1mw4/\nopn1lw1", "principal", "accessory", "full"))),  
  color = c("tan", "gold", "white", "violet", "deepskyblue", "chartreuse3", "orangered", "brown", "orange", "lightgrey"), 
  to.upper = FALSE)

# Species image files
# image.files = list.files("../../figures/animals/black", full.names = TRUE) %>% setNames(c('Chicken', 'Human', 'Lizard', 'Opossum', 'Rat', 'Zebrafish'))

major_annotation_pr_ze = CelltypeAnnotation(
  gene = list(c('isl2b', 'rbpms2a', 'rbpms2b', 'slc17a6a', 'slc17a6b', 'robo2'), 
              c('vsx1', 'vsx2', 'otx2', 'isl1', "prkcab", "prkcaa", 'grik1a', 'grik1b'),
              c('pax6a', 'tfap2a', 'tfap2b', 'tfap2c', 'gad1a', 'gad1b', 'gad2', 'slc6a9', 'slc6a5'), 
              c("onecut2"), 
              c('rho', 'rhol', 'opn1lw1','opn1mw1', 'opn1sw1','opn1sw2', 'arr3a', 'pde6h'),
              c("slc1a2b", "rlbp1b",  "apoeb", "apoba")), 
  annotation = factor(c("RGC",  # blue
                        "BC",  # chartreuse2
                        "AC", # cyan
                        "HC", # gold
                        "PR", # red
                        "MG"), # magenta
                      levels = c("PR", "HC", "BC", "AC", "RGC", "MG")),  
  color = c("blue", 
            "chartreuse2", 
            "cyan", 
            "gold", 
            "red", 
            "magenta"), to.upper = FALSE)

major_annotation_pr_ze_gaba_gly = CelltypeAnnotation(
  gene = list(c('isl2b', 'rbpms2a', 'rbpms2b', 'slc17a6a', 'slc17a6b', 'pou4f3', 'robo2'), 
              c('vsx1', 'vsx2', 'otx2', 'isl1', "prkcab", "prkcaa", 'grik1a', 'grik1b'),
              c('pax6a', 'tfap2a', 'tfap2b', 'tfap2c', 'gad2', 'gad1a', 'gad1b'), 
              c('slc6a9', 'slc6a5'), 
              c("onecut2"), 
              c('rho', 'rhol', 'opn1lw1','opn1mw1', 'opn1sw1','opn1sw2', 'arr3a', 'pde6h'),
              c("slc1a2b", "rlbp1b",  "apoeb", "apoba")), 
  annotation = factor(c("RGC",  # blue
                        "BC",  # chartreuse2
                        "gabaAC", # cyan
                        'glyAC',
                        "HC", # gold
                        "PR", # red
                        "MG"), # magenta
                      levels = c("PR", "HC", "BC", 'glyAC', "gabaAC", "RGC", "MG")),  
  color = c("blue", 
            "chartreuse2", 
            "cyan", 
            'cyan4',
            "gold", 
            "red", 
            "magenta"), to.upper = FALSE)

# pr_annotation = CelltypeAnnotation(
#   gene = list(c("PTPRD", "OPN1SW", "NRXN3", "RS1", "GNGT2"), 
#               c("THRB", "SOX5", "ACSL1", "OPN1MW", "OPN1LW")
#               c("RHO", "ESPRRB", "GNAT1", "PDE6B", "PDE6G", "PDC")),
#   annotation = factor(c("UV", 
#                         "Blue", 
#                         "Green", 
#                         ), # magenta
#                       levels = c("Rod", "ML_cone", "S_cone")),  
#   color = c("violet", "blue", "green", "red", "grey"))

# Consistent colors across plots; don't use globally just yet
# options(ggplot2.discrete.colour= list(dtColors()))

major_annotation_ml = CelltypeAnnotation(
  gene = list(c("RBPMS", "SLC17A6", "NEFL"), #"POU6F2", "RBFOX3"), 
              c('OTX2', 'GRM6', 'GRIK1'),  
              c("PAX6", "GAD1", "GAD2", "SLC6A9"), 
              c("ONECUT1", "LHX1", "CALB1"), 
              c('CRX', 'PDE6H', 'ARR3'), 
              c("SAG", "PDC", "RHO"), 
              c("SLC1A3", "RLBP1",  "APOE")), 
  annotation = factor(c("RGC",  # blue
                        "BC",  # chartreuse2
                        "AC", # cyan
                        "HC", # gold
                        "Cone",# orange
                        "Rod",  # red
                        "MG"), # magenta
                      levels = c("Rod", "Cone", "HC", "BC", "AC", "RGC", "MG")),  
  color = colorspace::lighten(c("blue", 
                                "chartreuse2", 
                                "cyan", 
                                "gold", 
                                "red", 
                                "grey", 
                                # "red",
                                "magenta"), 0.4))

major_annotation_rat = CelltypeAnnotation(
  gene = list(c("RBPMS", "SLC17A6", "NEFL"), #"POU6F2", "RBFOX3"), 
              c('OTX2', 'GRM6', 'GRIK1'),  
              c("PAX6", "GAD1", "GAD2", "SLC6A9"), 
              c("ONECUT1", "LHX1", "CALB1"), 
              c("ARR3", 'OPN1MW', 'OPN1SW'), 
              c("SAG", "PDC", "RHO"), 
              c("SLC1A3", "RLBP1",  "APOE")), 
  annotation = factor(c("RGC",  # blue
                        "BC",  # chartreuse2
                        "AC", # cyan
                        "HC", # gold
                        "Cone", # orange
                        "Rod",  # red
                        "MG"), # magenta
                      levels = c("Rod", "Cone", "HC", "BC", "AC", "RGC", "MG")),  
  color = colorspace::lighten(c("blue", 
                                "chartreuse2", 
                                "cyan", 
                                "gold", 
                                "red", 
                                "grey", 
                                # "red",
                                "magenta"), 0.4))


# 'NPY', 
# 'TFAP2D', 
# 'OPN4', 
# 'NOS1', 
# 'LOC116943951' # EOMES-like gene
major_annotation_lamprey = CelltypeAnnotation(
  gene = list(c("RBPMS", 'RBPMS2', "SLC17A6", 'POU4F3', "NEFL", "NEFM", 'TFAP2D', 'OPN4', 'LOC116943951'), #"POU6F2", "RBFOX3"), 
              c('VSX2', 'OTX2', 'ISL1', 'GRM6', 'GRIK1'), #"PRKCA"), 
              c("PAX6", "TFAP2A", "TFAP2B", "TFAP2C", "GAD1", "GAD2", 'NPY', 'NOS1'), 
              c("SLC6A9", 'SLC6A5'),
              c("ONECUT1", 'ONECUT2', 'ONECUT3', "CALB1", 'CALB2', "TPM3"), 
              # c("PDE6H", "CRX", "ARR3"), 
              # c("SAG", "PDC", "RHO"), 
              c("PDE6H", "CRX", "ARR3","SAG", "PDC", "RHO"),
              c("SLC1A3", "RLBP1",  "APOE", "APOB")), 
  annotation = factor(c("RGC",  # blue
                        "BC",  # chartreuse2
                        "gabaAC", # cyan
                        'glyAC',
                        "HC", # gold
                        # "Cone",# orange
                        # "Rod",  # red
                        "PR", # red
                        "MG"), # magenta
                      levels = c("PR", "BC", "HC", 'glyAC', "gabaAC", "RGC", "MG")),  
  color = c("blue", 
            "chartreuse2", 
            "cyan", 
            'cyan4',
            "gold", 
            "red", 
            # "darkred", 
            # "red",
            "magenta"))

major_annotation_lamprey_on_off = CelltypeAnnotation(
  gene = list(c("RBPMS", 'RBPMS2', "SLC17A6", 'POU4F3', "NEFL", "NEFM", 'TFAP2D', 'OPN4', 'LOC116943951'), #"POU6F2", "RBFOX3"),
              c('VSX2', 'OTX2', 'ISL1', 'GRM6'),
              c('GRIK1'),
              c("PAX6", "TFAP2A", "TFAP2B", "TFAP2C", "GAD1", "GAD2", 'NPY', 'NOS1'),
              c("SLC6A9", 'SLC6A5'),
              c("ONECUT1", 'ONECUT2', 'ONECUT3', "CALB1", 'CALB2', "TPM3"),
              # c("PDE6H", "CRX", "ARR3"),
              # c("SAG", "PDC", "RHO"),
              c("PDE6H", "CRX", "ARR3","SAG", "PDC", "RHO"),
              c("SLC1A3", "RLBP1",  "APOE", "APOB")),
  annotation = factor(c("RGC",  # blue
                        "onBC",  # chartreuse2
                        "offBC",
                        "gabaAC", # cyan
                        'glyAC',
                        "HC", # gold
                        # "Cone",# orange
                        # "Rod",  # red
                        "PR", # red
                        "MG"), # magenta
                      levels = c("PR", "onBC", "offBC", "HC", 'glyAC', "gabaAC", "RGC", "MG")),
  color = c("blue",
            "chartreuse1", 
            'chartreuse4',
            "cyan",
            'cyan4',
            "gold",
            "red",
            # "darkred",
            # "red",
            "magenta"))

major_annotation_pr_gaba_gly = CelltypeAnnotation(
  gene = list(c("RBPMS", 'RBPMS2', "SLC17A6", 'POU4F3', "NEFL", "NEFM"),
              c('VSX2', 'VSX1', 'OTX2', 'ISL1', 'GRM6', 'GRIK1', "PRKCA"), 
              c("PAX6", "TFAP2A", "TFAP2B", "TFAP2C", "GAD1", "GAD2"),
              c("SLC6A9", 'SLC6A5'), 
              c("ONECUT1", "LHX1", "CALB1", "TPM3"), 
              c("PDE6H", "CRX", "ARR3","SAG", "PDC", "RHO"),
              c("SLC1A3", "RLBP1",  "APOE", "APOB")), 
  annotation = factor(c("RGC",  # blue
                        "BC",  # chartreuse2
                        "gabaAC", # cyan
                        'glyAC',
                        "HC", # gold
                        "PR", # red
                        "MG"), # magenta
                      levels = c("PR", "BC", "HC", 'glyAC', "gabaAC", "RGC", "MG")),  
  color = c("blue", 
            "chartreuse2", 
            "cyan", 
            'cyan4',
            "gold", 
            "red",
            "magenta"))

major_annotation_pr_gaba_gly_on_off = CelltypeAnnotation(
  gene = list(c("RBPMS", 'RBPMS2', "SLC17A6", 'POU4F3', "NEFL", "NEFM"),
              c('VSX2', 'OTX2', 'ISL1', 'GRM6'),
              c('GRIK1'), 
              c("PAX6", "TFAP2A", "TFAP2B", "TFAP2C", "GAD1", "GAD2"),
              c("SLC6A9", 'SLC6A5'), 
              c("ONECUT1", "LHX1", "CALB1", "TPM3"), 
              c("PDE6H", "CRX", "ARR3","SAG", "PDC", "RHO"),
              c("SLC1A3", "RLBP1",  "APOE", "APOB")), 
  annotation = factor(c("RGC",  # blue
                        "onBC",  # chartreuse2
                        'offBC',
                        "gabaAC", # cyan
                        'glyAC',
                        "HC", # gold
                        "PR", # red
                        "MG"), # magenta
                      levels = c("PR", "HC", "onBC", 'offBC', 'glyAC', "gabaAC", "RGC", "MG")),  
  color = c("blue", 
            "chartreuse1", 
            'chartreuse4',
            "cyan", 
            'cyan4',
            "gold", 
            "red",
            "magenta"))

# GetAnnotationColors = function(object, annotation = major_annotation_pr_gaba_gly, group.by = 'annotated', color.by = 'cell_class2'){
#   metadata = Metadata(object, group.by, color.by)
#   metadata$colors = ColorCode(annotation, metadata[[color.by]])
#   return(metadata$colors %>% setNames(metadata[[group.by]]))
# }

# more palettes
major_annotation_palette = c(unique(Colors(major_annotation_pr)) %>% setNames(unique(Annotation(major_annotation_pr))), 'Other' = 'grey')
major_annotation_palette2 = c(unique(Colors(major_annotation_pr_gaba_gly)) %>% setNames(unique(Annotation(major_annotation_pr_gaba_gly))), 'Other' = 'grey')
# major_annotation_palette3 = c(major_annotation_palette2, AC = 'cyan')
major_annotation_palette3 = c(major_annotation_palette2, AC = 'cyan', Rod = 'darkred', Cone = 'red')

cell_class2_colors = ColorCode(major_annotation_pr_gaba_gly, 
                               unique(Annotation(major_annotation_pr_gaba_gly))) %>%
  setNames(unique(Annotation(major_annotation_pr_gaba_gly)))

cell_class3_colors = ColorCode(major_annotation_pr_gaba_gly_on_off, 
                               unique(Annotation(major_annotation_pr_gaba_gly_on_off))) %>%
  setNames(unique(Annotation(major_annotation_pr_gaba_gly_on_off)))

index = read.csv('../AC/samap_index.csv')
CELL_CLASSES = setdiff(names(major_annotation_palette), 'Other') %>%
  setNames(setdiff(names(major_annotation_palette), 'Other'))
CELL_CLASSES2 = setdiff(names(major_annotation_palette2), 'Other') %>%
  setNames(setdiff(names(major_annotation_palette2), 'Other'))
CELL_CLASSES3 = c(names(cell_class3_colors), 'Other') %>% setNames(c(names(cell_class3_colors), 'Other'))

NM_SPECIES = c('Chicken', 'Lizard', 'Axolotl', 'Newt', 'Zebrafish', 'Killifish', 'Goldfish', 'Shark', 'Lamprey') %>% 
  setNames(c('Chicken', 'Lizard', 'Axolotl', 'Newt', 'Zebrafish', 'Killifish', 'Goldfish', 'Shark', 'Lamprey'))
nm_palette = species_palette3[NM_SPECIES]
nm_palette2 = nm_palette %>% setNames(index$ident[match(names(nm_palette), index$species)])

# Type colors
if(!exists('ac.ortho')) ac.ortho = LoadACOrtho()
dend_custom_order = readRDS('../../Ortho_Objects/ac.dendrogram.v2.rds')
clusters <- cutree(as.hclust(dend_custom_order), k = 10)
type_cols2.no = colorspace::lighten(ClusterPalette(sort(clusters), 
  colors = c('deeppink', 'cyan', 'orange', 'gold', 'chartreuse', '#584B9FFF', 'blue', 'tan', 'deepskyblue', 'violet'), vary.by = 0.4), 0.2)
names(type_cols2.no) = as.numeric(names(sort(clusters))) + 1
type_cols2.no = type_cols2.no[order(as.numeric(names(type_cols2.no)))]
type_cols2 = type_cols2.no %>% setNames(levels(ac.ortho$type_unordered))
type_cols3 = type_cols2 # New orthotype palette
names(type_cols3) = orthotype_labels(names(type_cols2))
type_cols3.no = type_cols3
names(type_cols3.no) =  gsub('oAC|\\*', '', ExtractString(names(type_cols3), after = ' '))

# Wider gamut palette
dend_custom_order = readRDS('../../Ortho_Objects/ac.dendrogram.v2.rds')
labels(dend_custom_order) = orthotype_no_labels(as.numeric(labels(dend_custom_order)) + 1)
clusters <- cutree(as.hclust(dend_custom_order), k = 13)
clusters = clusters[as.character(labels(dend_custom_order))]
type_cols4.no = ClusterPalette2(clusters,
                               colors = c('blue', '#584B9FFF',   'deepskyblue', 'cyan', 'brown','tan', 'red', 
                                          'orange', 'deeppink', 'chartreuse', 'antiquewhite', 'violet', 'gold'), 
                               debug = F)
names(type_cols4.no) = names(clusters) # orthotype_no_labels(as.numeric(names(sort(clusters)))+1)
type_cols4 = type_cols4.no
names(type_cols4) = levels(ac.ortho$orthotype)

# Paths
hahn.path = '~/Google Drive/Shared drives/Shekharlab_data3/Josh H Projects/Evolution of Cell Types/Analysis/Species_Objects/Final Objects/'

# General AC objects
neuropeptides = read.table("../../reference_files/complete_neuropeptide_symbols.txt")$V1
hk_genes = read.table("../../reference_files/Eisenberg_Cell_2013_hk_genes.txt")$V1
tf_genes = read.csv("../../reference_files/cisbp_20220920.csv")$Name

# Global variables
NCLUSTERS = 42
type_cols = AcPalette2(42) %>% setNames(1:NCLUSTERS)
speciesList = names(SPECIESFILES)
speciesList.no.gf = setdiff(names(SPECIESFILES), 'Goldfish')
speciesList.liz = speciesList[1:17] #setdiff(speciesList.no.gf, c('Zebrafish', 'Lamprey')) %>% setNames(setdiff(speciesList.no.gf, c('Zebrafish', 'Lamprey')))
speciesList.mammals = speciesList.liz[1:15]
JACCARD.CUTOFF = 0.25
KNOWNTYPES = c('SAC', 'A2', 'VG3') %>% setNames(c('SAC', 'A2', 'VG3'))

# Global files
neuropeptides = read.table("../../reference_files/complete_neuropeptide_symbols.txt")$V1
hk_genes = read.table("../../reference_files/Eisenberg_Cell_2013_hk_genes.txt")$V1
tf_genes = read.csv("../../reference_files/cisbp_20220920.csv")$Name

# Feature plot colors
fcols = c('grey', 'deeppink')

NCLUSTERS = 42
type_cols = AcPalette2(42) %>% setNames(1:NCLUSTERS)

SPECIES_PALETTE = myPalette(length(SPECIESFILES)) %>% setNames(names(SPECIESFILES))
speciesList.cp = c('Lamprey', 'Shark', 'Goldfish', 'Killifish', 'ZebrafishLyu', #'ZebrafishLyuPR', 
                   'Axolotl', 'Newt', 
                   'Chicken', 'Lizard', 'Opossum', 'Squirrel', 'Rat', 
                   'TreeShrew', #'TreeShrewClean', 
                   'MouseLemur', 'Marmoset', 'Macaque', 'Human')

species_metadata <- data.frame(
  Species = c(
    "Human", "Macaque", "Marmoset", "MouseLemur", "TreeShrew",
    "Mouse", "Rhabdomys", "Rat", "Peromyscus", "Squirrel",
    "Cow", "Sheep", "Pig", "Ferret", "Opossum", "Chicken", "Lizard"
  ),
  Activity = c(
    "Diurnal", "Diurnal", "Diurnal", "Nocturnal", "Diurnal",
    "Nocturnal", "Diurnal", "Nocturnal", "Nocturnal", "Diurnal",
    "Diurnal", "Diurnal", "Diurnal", NA, "Nocturnal", "Diurnal", "Diurnal"
  ), 
  Clade = c('Primate', 'Primate', 'Primate', 'Primate', 'Scandentia', 
            'Rodent', 'Rodent', 'Rodent', 'Rodent', 'Rodent', 
            'Laurasiatheria', 'Laurasiatheria', 'Laurasiatheria', 'Laurasiatheria', 
            'Marsupial', 'Sauropsid', 'Sauropsid')
  )

species_seq_method <- data.frame(
  CommonName = c(
    "Lamprey", "Zebrafish", "Goldfish", "Lizard", "Chicken",
    "Opossum", "Ferret", "Pig", "Cow", "Sheep",
    "Squirrel", "Peromyscus", "Rat",
    "Rhabdomys", "Mouse", "TreeShrew",
    "MouseLemur", "Marmoset", "Macaque", "Human"
  ),
  Method = c(
    "scRNA-seq", "snRNA-seq", "scRNA-seq", "snRNA-seq", "scRNA-seq",
    "snRNA-seq", "snRNA-seq", "both", "snRNA-seq", "snRNA-seq",
    "snRNA-seq", "scRNA-seq", "snRNA-seq", "snRNA-seq", "scRNA-seq", "snRNA-seq",
    "snRNA-seq", "scRNA-seq", "scRNA-seq", "snRNA-seq"
  ),
  stringsAsFactors = FALSE
)


ANTIBODIES <- unique(c(
  "CACNA1S","CACNA1S","ACHE","ACHE","ACT","ACTB","ACTB","ACTB",
  "ALDH1A2","ALDH1A2","AMIGO1","AMIGO1","AMIGO1","ANK2","ANK2",
  "ANK3","ANK3","ANK3","ATP5F1A","BNC2","BNC2","BNC2","BNC2",
  "BNC2","BNC2","BNC2","BNC2","BNC2","BNC2","BNC2","BNC2",
  "POU4F1","POU4F1","POU4F1","POU4F2","POU4F2","POU4F2",
  "POU4F3","POU4F3","CA8","CA8","CA8","CA8","CA10","CA10",
  "CABP1","CABP5","CABP5","CABP5","CABP5","CDH13","CALB1",
  "CALB1","CALB1","CALB1","CALB1","CALB1","CALB1","CALB1",
  "CALB1","CALB1","CALB1","CALB2","CALB2","CALB2","CALB2",
  "CAMK2A","CA8","CA8","CARTPT","CARTPT","CARTPT","CARTPT",
  "CACNA1C","CACNA1C","CACNA1D","CACNA1D","CACNA1G",
  "CACNA1G","CACNA1H","CACNA1H","CBFB","CBFB","CBFB",
  "CBLN1","CBLN1","MCAM","MCAM","MCAM","FUT4","PECAM1",
  "CD68","CD68","CD83","CBLN2","CBLN2","CALCA","CALCA",
  "CALCA","CALCA, CALCB","CALCA, CALCB","CALCA","DLG2",
  "DLG2","DLG2","DLG2","DLG2","CHAT","CHAT","CHAT","CHAT",
  "CCK","CCK","CHRNA6","CBX1","VSX2","SHISA9","SHISA9",
  "CNGB1","CNGB1","CNGB1","COCH","COCH","COCH","ARR3",
  "GJC1","CNIH1","Not found","CALCRL","CTBP","CtBP",
  "CYP26A1","DRD1","DAB1","DAB1","DAB1","DSCAM","SLC1A2",
  "EFHD2","ELFN2","ELFN2","EOMES","RAPGEF4","RAPGEF4",
  "RAPGEF4","RAPGEF4","EPS8","ESRRA","EXT1","EXT1","FEZF1",
  "FGF14","FGF14","FGF7","FGF7","FOXP1","FOXP1","FOXP1",
  "FOXP2","FOXP2","FOXP2","FOXP2","FOXP3","FOXP3","FOXP4",
  "FSTL4","FSTL4","FXR2","GABRB2","GABRB2","GABRB2",
  "GABRG3","GABRG3","GAD1 / GAD2","GAD1 / GAD2","GAD1",
  "GNAO1","GNAO1","GAPDH","GBX2","GBX2","GBX2","GBX2",
  "GBX2","GPHN","GPHN","GPHN","GFAP","GFAP","Not found",
  "SLC1A2","GRIA1","GRIA1","GRIA1","GRIA1","GRIA1",
  "GRIA2","GRIA2","GRIA2","GRIA2","GRIA2","GRIA2",
  "GRIA2","GRIA3","GRIA3","GRIA4","GRIA4","GRIA4",
  "GRIK1","GRIK1","GRIK1","GRIK1","GRIK1","GRIK1",
  "GRIK1","GRIK1/2/3","GRIK1/2/3","GRIK1/2/3","GRIK2",
  "GRIK2","GRIK2","GRIK2/3","GRIK2","GLUL","GLUL",
  "GLRA1","GLRA1","GLRA3","GLRA3","Alias symbol GLYT1",
  "Alias symbol GLYT1","SLC6A9","GNAT2","GNAT2","GNAT2",
  "GNG13","GNG13","GRIK1","GRIK1","GRIP1","GRIP1","HA",
  "CD44","HCN1","HCN1","HCN2","HCN2","HCN3","HCN3",
  "HCN4","HCN4","HAL","AIF1","AIF1","AIF1","IGFBP5",
  "ITPR1","ITPR1","ITPR1","ITPR1","ISL1 & ISL2",
  "ISL1 & ISL2","ISL2","GRIK4","GRIK5","GRIK5","GRIK5",
  "KCNIP1","MKI67","KIF2A","KIF2A","KCNJ2","KCNJ2",
  "KCNJ10","KCNJ10","XRCC5","KCNA1","KCNA1","KCNA1",
  "KCNA1","KCNA1","KCNA1","KCNA2","KCNA2","KCNA3",
  "KCNA3","KCNA4","KCNA4","KCNA5","KCNA5","KCNA6",
  "KCNA6","KCNB1","KCNB1","KCNB1","KCNB1","KCNB1",
  "KCNB1","KCNB2","KCNB2","KCNB2","KCNB2","KCNB2",
  "KCNB2","KCNB2","KCNB2","KCNB2","KCNB2","KCNB2",
  "KCNC1","KCNC1","KCNC3","KCNC3","HCNC4","HCNC4",
  "KCND2","KCND2","KCND3","KCND3","KCNG4","KCNG4",
  "KCNQ1","KCNQ2","KCNQ4","KCNV2","KCNV2","KCNV2",
  "KCNV2","KCNV2","KCNV2","KCNV2","KCNV2","KCNV2",
  "LRIT3","LRIT3","LRIT3","LRIT3","LYPD1","MAF","MAF",
  "MAF","MAF","MAFB","MAFB","MARCKS","MEIS2","MEIS2",
  "MEIS2","MEIS2","OPN4","OPN4","OPN4","OPN4","OPN4",
  "GRM4","GRM6","GRM6","GRM6","GRM6","GRM6","GRM6",
  "GRM6","GRM6","GRM6","GRM7","TRPM1","HSPA9","MYLK",
  "MYLK","ATP1A3","ATP1A3","SCN1A","SCN1A","SCN2A",
  "SCN2A","SCN8A","SCN8A","SCN4B","SCN4B","CDH2",
  "CDH2","NETO1","NETO1","NETO1","NETO2","NETO2",
  "RBFOX3","RBFOX3","NEUROD2","NEUROD2","NEUROD2",
  "NEFH","NEFH","NEFM","NPY","NPY","NTS","NTS","NTS",
  "TACR3","TACR3","TACR3","TACR3","NOS1","NOS1",
  "NR2F2","NR2F2","NXPH1","NXPH1","NXPH2","ONECUT1",
  "ONECUT1","ONECUT1","OPN1LW / OPN1MW","OPN1SW",
  "P2RX4","P2RX6","KCNIP1/2/3/4","KCNIP1/2/3/4",
  "KCNAB1/2","KCNAB1/2","NFASC","NFASC","PVALB",
  "PVALB","PVALB","PVALB","PVALB","PAX6","PAX6",
  "PAX6","PAX6","PCP2","PCP2","PDE4B","PDE4B",
  "PDGFRA","PDGFRA","PDGFRA","PDGFRA","PDGFRA",
  "PDGFRA","PDGFRA","PIKFYVE","PRKAR2B","PRKAR2B",
  "PRKCA","PRKCA","PRKCA","PRKCA","PRKCA","PRKCA",
  "PLCB4","PLCB4","POU6F2","POU6F2","PROX1","PROX1",
  "PRPH2","PRPH2","DLG4","DLG4","DLG4","DLG4","DLG4",
  "QPCT","RARA","RARA","RARB","RARB","RARB","RARB",
  "RBPMS","RBPMS","RBPMS","RBPMS","RBPMS","RBPMS",
  "RCVRN","RCVRN","RCVRN","RCVRN","RCVRN","RCVRN",
  "RELN","RHO","RHO","CTBP2","CTBP2","CTBP2",
  "RPE65","RPE65","RUNX1","RUNX1","RUNX1","RUNX1",
  "RUNX1","RYR1","RYR1","SATB2","SATB2","SATB2",
  "SATB2","SATB2","SATB2","SATB2","SCGN","SCGN",
  "HTOR","SGK1","SLC35D3","SLC35D3","SLC35D3",
  "SLC35D3","KCNMA1","SMYD2","SOX5","SOX5","SPON1",
  "Stxbp6","Stxbp6","SYP","SYT2","SYT10","SYT10",
  "SYT3","SYT3","SYT5","SYT5","SYT6","SYT6","SYT6",
  "SYT7","SYT7","SYT9","SYT9","SYT1","CACNG2",
  "CACNG2","CACNG2","CACNG2","CACNG2/4/8",
  "CACNG2/4/8","CACNG3","CACNG8","TBR1","TBR1",
  "EOMES","EOMES","TCF7L2","TCF7L2","TCF7L2",
  "TCF7L2","TCF7L2","TCF7L2","TCF7L2","TFAP2D",
  "TFAP2D","TGFBI","THY1","KCNK4","KCNK4","TRHDE",
  "TRPM1","TRPM1","TRPM1","TRPM1","TUBB","TRARG1",
  "TH","TH","VIP","VIP","SLC32A1","SLC17A7",
  "SLC17A7","SLC17A7","SLC17A7","SLC17A7",
  "SLC17A8","SLC17A8","SLC17A8","SLC17A8",
  "SLC17A8","CACNA1F","CACNA2D4","ZBTB16",
  "TJP1","TJP1"
))

