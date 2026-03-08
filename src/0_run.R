parameters = data.frame(species = phylogenetic_order, 
                        group.by = c("animal", "animal", "animal", "animal", "animal", "batch", "animal", "animal", "animal", "animal", "animal", "animal", "animal", "animal", "animal", "animal", "animal", "orig.file", "animal", "animal"), 
                        method = c("seurat", "seurat", "harmony", "harmony", "seurat", "harmony", "seurat", "harmony", "harmony", "harmony", "seurat", "seurat", "harmony", "seurat", "seurat", "harmony", "seurat", "harmony", "seurat", "harmony"))

# Class and type annotation for new channels
rmarkdown::render("ClassTypeAnnotation.Rmd", output_file = "html_reports/ClassTypeAnnotation.html", params = list(species = "Squirrel",
                                                                                                                  species_dir = "/clusterfs/kslab/DATA/BROAD_DATA/Squirrel2/",
                                                                                                                  samples = list("sq14percneun_s1", "sq14percneun_s2"), 
                                                                                                                  add = TRUE,
                                                                                                                  metadata = TRUE,
                                                                                                                  save = TRUE))

rmarkdown::render("ClassTypeAnnotation.Rmd", output_file = "html_reports/ClassTypeAnnotation.html", params = list(species = "Sheep",
                                                                                                                  species_dir = "/clusterfs/kslab/DATA/BROAD_DATA/Lamb/",
                                                                                                                  samples = list("sh07percneun_s1", "sh07percneun_s2"),
                                                                                                                  add = TRUE,
                                                                                                                  metadata = TRUE,
                                                                                                                  save = TRUE))

rmarkdown::render("ClassTypeAnnotation.Rmd", output_file = "html_reports/ClassTypeAnnotation.html", params = list(species = "MouseLemur",
                                                                                                                  species_dir = "/clusterfs/kslab/DATA/BROAD_DATA/Lamb/",
                                                                                                                  samples = list("sh07percneun_s1", "sh07percneun_s2"),
                                                                                                                  add = TRUE,
                                                                                                                  metadata = TRUE,
                                                                                                                  save = TRUE))

rmarkdown::render("ClassTypeAnnotation.Rmd", output_file = "html_reports/ClassTypeAnnotation.html", params = list(species = "Rat",
                                                                                                                  species_dir = "/clusterfs/kslab/DATA/BROAD_DATA/Lamb/",
                                                                                                                  samples = list("sh07percneun_s1", "sh07percneun_s2"),
                                                                                                                  add = TRUE,
                                                                                                                  metadata = TRUE,
                                                                                                                  save = TRUE))
# Redid 11/23
render_report(species = "Human",
              initial = TRUE,
              batch_int = TRUE,
              integrate_by = "animal",
              harmony = FALSE,
              contamination_threshold = 6, 
              contamination = list(rod = NA, cone = NA, hc = NA, bc = NA, ac = NA, rgc = NA, mg = NA), 
              nFeature_threshold = -2.5,
              doublet_finder = TRUE,
              manual_annotation = list(SAC = c(10), A17 = c(NA), nNOS = c(31,34), CA1 = c(23), CA2 = c(11), VG1 = c(NA), A2 = c(3), SEG = c(9), VG3 = 6), 
              sac_annotation = FALSE,
              de_expression = TRUE, 
              save = TRUE)

render_report(species = "Macaque",
              initial = TRUE,
              batch_int = TRUE,
              integrate_by = "animal",
              harmony = FALSE,
              contamination_threshold = 6, 
              contamination = list(rod = NA, cone = NA, hc = NA, bc = NA, ac = NA, rgc = NA, mg = NA), 
              nFeature_threshold = -2.5,
              doublet_finder = FALSE,
              manual_annotation = list(SAC = c(18), A17 = c(1,2,3), nNOS = c(35), CA1 = c(43), CA2 = c(11), VG1 = c(NA), A2 = c(15), SEG = c(0,8,16), VG3 = 34), 
              sac_annotation = FALSE,
              de_expression = TRUE,
              save = TRUE)

# Checked 8/14/23
render_report(species = "Marmoset",
              initial = TRUE,
              batch_int = TRUE,
              integrate_by = "animal",
              harmony = TRUE,
              contamination_threshold = 6, 
              contamination = list(rod = c(NA), cone = c(NA), hc = c(NA), bc = c(13,26), ac = NA, rgc = c(NA), mg = c(NA)), 
              nFeature_threshold = -2.5,
              doublet_finder = TRUE,
              manual_annotation = list(SAC = c(15), A17 = c(NA), nNOS = c(24), CA1 = c(21), CA2 = c(0), VG1 = c(NA), A2 = c(5), SEG = c(6,7,17), VG3 = 14), 
              sac_annotation = FALSE,
              de_expression = TRUE, 
              save = TRUE)

render_report(species = "MouseLemur",
              initial = TRUE,
              batch_int = FALSE,
              integrate_by = "animal",
              harmony = FALSE,
              contamination_threshold = 6, 
              contamination = list(rod = c(NA), cone = c(NA), hc = c(NA), bc = c(20,24,43,44,45), ac = NA, rgc = NA, mg = c(NA)), 
              nFeature_threshold = -2.5,
              doublet_finder = FALSE,
              manual_annotation = list(SAC = c(33), A17 = NA, nNOS = c(17,19), CA1 = c(36), CA2 = c(NA), VG1 = c(NA), A2 = c(3), SEG = c(4,5,8), VG3 = 0), 
              sac_annotation = FALSE,
              de_expression = TRUE, 
              save = FALSE)

render_report(species = "TreeShrew",
              initial = TRUE,
              batch_int = TRUE,
              integrate_by = "animal",
              harmony = TRUE,
              contamination_threshold = 6, 
              contamination = list(rod = c(NA), cone = c(7), hc = c(NA), bc = c(NA), ac = c(37), rgc = c(56,70), mg = c(NA)), 
                              # list(rod = c(NA), cone = c(8,16,45), hc = c(NA), bc = c(NA), ac = NA, rgc = c(55,66), mg = c(38))
              nFeature_threshold = -2,
              doublet_finder = TRUE,
              manual_annotation = list(SAC = c(0), A17 = c(NA), nNOS = c(5,7), CA1 = c(50), CA2 = c(NA), VG1 = c(NA), A2 = c(48), SEG = c(18), VG3 = 14), 
              sac_annotation = FALSE,
              de_expression = TRUE,
              save = TRUE)

render_report(species = "Rhabdomys",
              initial = TRUE,
              batch_int = FALSE,
              integrate_by = "animal",
              harmony = FALSE,
              contamination_threshold = 6, 
              contamination = list(rod = c(2,32,33,36), cone = c(NA), hc = c(NA), bc = c(9,34,38), ac = NA, rgc = c(NA), mg = c(NA)),
              nFeature_threshold = -2.5,
              doublet_finder = FALSE,
              manual_annotation = list(SAC = c(0,9), A17 = c(NA), nNOS = c(NA), CA1 = c(24), CA2 = c(23), VG1 = c(NA), A2 = c(6), SEG = c(1), VG3 = 2), 
              sac_annotation = FALSE,
              de_expression = TRUE, 
              save = TRUE)

render_report(species = "Rat",
              initial = TRUE,
              batch_int = FALSE,
              integrate_by = "animal",
              harmony = TRUE,
              contamination_threshold = 6, 
              contamination = list(rod = c(NA), cone = c(NA), hc = c(NA), bc = c(40,46,58,60,61), ac = NA, rgc = NA, mg = NA), 
              nFeature_threshold = -3,
              doublet_finder = FALSE,
              manual_annotation = list(SAC = c(10,14), A17 = 1, nNOS = c(20,52,53), CA1 = c(54), CA2 = c(6,22,49), VG1 = c(NA), A2 = c(3), SEG = c(4,11), VG3 = 35), 
              sac_annotation = FALSE,
              de_expression = TRUE, 
              save = FALSE)

render_report(species = "Peromyscus",
              initial = TRUE,
              batch_int = TRUE,
              integrate_by = "animal",
              harmony = TRUE,
              contamination_threshold = 6, 
              contamination = list(rod = c(NA), cone = c(NA), hc = c(NA), bc = c(5,17,19), ac = NA, rgc = c(18), mg = c(NA)),
              nFeature_threshold = -2.5,
              doublet_finder = TRUE,
              manual_annotation = list(SAC = c(2), A17 = c(3), nNOS = c(5,32,36), CA1 = c(33), CA2 = c(15), VG1 = c(NA), A2 = c(1), SEG = c(0,7,34), VG3 = 28),
              sac_annotation = FALSE,
              de_expression = TRUE, 
              save = TRUE)

render_report(species = "Squirrel",
              initial = TRUE,
              batch_int = TRUE,
              integrate_by = "animal",
              harmony = TRUE,
              contamination_threshold = 6, 
              contamination = list(rod = c(NA), cone = c(NA), hc = c(NA), bc = c(23,43,48,52), ac = NA, rgc = c(NA), mg = c(45)), 
              nFeature_threshold = -1.85,
              doublet_finder = TRUE,
              manual_annotation = list(SAC = c(0,1), A17 = c(NA), nNOS = c(13,18,21,23,37), CA1 = c(NA), CA2 = c(40), VG1 = c(NA), A2 = c(45), SEG = c(7), VG3 = 2), 
              sac_annotation = FALSE,
              de_expression = TRUE, 
              save = TRUE)

render_report(species = "Ferret",
              initial = TRUE,
              batch_int = TRUE,
              integrate_by = "animal",
              harmony = FALSE,
              contamination_threshold = 4, 
              contamination = list(rod = c(8), cone = c(NA), hc = c(NA), bc = c(21), ac = NA, rgc = c(NA), mg = c(NA)),
              nFeature_threshold = -2,
              doublet_finder = TRUE,
              manual_annotation = list(SAC = c(NA), A17 = c(NA), nNOS = c(6), CA1 = c(17), CA2 = c(7), VG1 = c(NA), A2 = c(16), SEG = c(8,14), VG3 = 13),
              sac_annotation = FALSE,
              de_expression = TRUE, 
              save = TRUE)

render_report(species = "Pig",
              initial = TRUE,
              batch_int = FALSE,
              integrate_by = "animal",
              harmony = FALSE,
              contamination_threshold = 4, 
              contamination = list(rod = c(13), cone = c(NA), hc = c(NA), bc = c(16), ac = NA, rgc = c(NA), mg = c(NA)), 
              nFeature_threshold = -2.5,
              doublet_finder = FALSE,
              manual_annotation = list(SAC = c(13), A17 = c(NA), nNOS = c(2,12), CA1 = c(21), CA2 = c(NA), VG1 = c(19), A2 = c(0), SEG = c(7,25), VG3 = 27),
              sac_annotation = FALSE,
              de_expression = TRUE, 
              save = FALSE)

render_report(species = "Cow",
              initial = TRUE,
              batch_int = FALSE,
              integrate_by = "animal",
              harmony = TRUE,
              contamination_threshold = 4, 
              contamination = list(rod = c(23,52,54), cone = c(NA), hc = c(NA), bc = c(36), ac = NA, rgc = c(NA), mg = c(NA)), 
              nFeature_threshold = -2.5,
              doublet_finder = FALSE,
              manual_annotation = list(SAC = c(1), A17 = c(NA), nNOS = c(17,19,34), CA1 = c(43), CA2 = c(36), VG1 = c(NA), A2 = c(2), SEG = c(15,3), VG3 = 33),
              sac_annotation = FALSE,
              de_expression = TRUE,
              save = TRUE)

render_report(species = "Sheep",
              initial = TRUE,
              batch_int = TRUE,
              integrate_by = "animal",
              harmony = FALSE,
              contamination_threshold = 6, 
              contamination = list(rod = c(43,49), cone = c(NA), hc = c(NA), bc = c(NA), ac = NA, rgc = c(NA), mg = c(NA)),
              nFeature_threshold = -2.5,
              doublet_finder = TRUE,
              manual_annotation = list(SAC = c(6), A17 = c(NA), nNOS = c(19,33), CA1 = c(38), CA2 = c(NA), VG1 = c(NA), A2 = c(20), SEG = c(23,28), VG3 = 43),
              sac_annotation = FALSE,
              de_expression = TRUE, 
              save = TRUE)

render_report(species = "Opossum",
              initial = TRUE,
              batch_int = TRUE,
              integrate_by = "animal",
              harmony = FALSE,
              contamination_threshold = 6, 
              contamination = list(rod = c(0,17,59), cone = c(NA), hc = c(NA), bc = c(NA), ac = NA, rgc = c(4,7), mg = c(NA)), 
              nFeature_threshold = -2.5,
              doublet_finder = TRUE,
              manual_annotation = list(SAC = c(1,3,46), A17 = NA, nNOS = c(44), CA1 = c(55), CA2 = c(NA), VG1 = c(NA), A2 = c(0), SEG = c(8,19,50), VG3 = 9), 
              sac_annotation = FALSE,
              de_expression = TRUE,
              save = TRUE)

render_report(species = "Chicken",
              initial = TRUE,
              batch_int = TRUE,
              integrate_by = "animal", # each channel was separate biological replicate
              harmony = TRUE,
              contamination_threshold = 6, 
              contamination = list(rod = NA, cone = NA, hc = NA, bc = 16, ac = NA, rgc = NA, mg = NA), 
              nFeature_threshold = -2.5,
              doublet_finder = TRUE,
              manual_annotation = list(SAC = 3, A17 = NA, nNOS = 54, CA1 = 6, CA2 = 5, VG1 = NA, A2 = 1, SEG = c(10,28,52), VG3 = 33), 
              sac_annotation = FALSE,
              de_expression = TRUE, 
              save = TRUE)

render_report(species = "Chicken_reclustered",
              initial = TRUE,
              batch_int = TRUE,
              integrate_by = "animal", # each channel was separate biological replicate
              harmony = TRUE,
              contamination_threshold = 6, 
              contamination = list(rod = NA, cone = c(38), hc = c(NA), bc = c(13,18), ac = NA, rgc = NA, mg = NA), 
              nFeature_threshold = -2.5,
              doublet_finder = TRUE,
              manual_annotation = list(SAC = c(2), A17 = c(NA), nNOS = c(63), CA1 = c(17), CA2 = c(NA), VG1 = c(NA), A2 = c(1), SEG = c(14,22), VG3 = 40), 
              sac_annotation = FALSE,
              de_expression = TRUE, 
              save = TRUE)

render_report(species = "Lizard",
              initial = TRUE,
              batch_int = TRUE,
              integrate_by = "animal",
              harmony = FALSE,
              contamination_threshold = 6, 
              contamination = list(rod = c(NA), cone = c(22), hc = c(NA), bc = c(33), ac = NA, rgc = c(NA), mg = c(NA)), 
              nFeature_threshold = -3,
              doublet_finder = TRUE,
              manual_annotation = list(SAC = c(12,24), A17 = c(NA), nNOS = c(35), CA1 = c(33), CA2 = c(NA), VG1 = c(NA), A2 = c(15), SEG = c(10,11), VG3 = 21), 
              sac_annotation = FALSE,
              de_expression = TRUE, 
              save = TRUE)

render_report(species = "Lizard_ncbi",
              initial = TRUE,
              batch_int = TRUE,
              integrate_by = "animal",
              harmony = FALSE,
              contamination_threshold = 6, 
              contamination = list(rod = c(NA), cone = c(22), hc = c(NA), bc = c(NA), ac = NA, rgc = c(NA), mg = c(NA)), 
              nFeature_threshold = -3,
              doublet_finder = TRUE,
              manual_annotation = list(SAC = c(14, 24), A17 = c(NA), nNOS = c(22, 18), CA1 = c(25), CA2 = c(NA), VG1 = c(NA), A2 = c(17), SEG = c(8,15), VG3 = 21), 
              sac_annotation = FALSE,
              de_expression = TRUE, 
              save = TRUE)

render_report(species = "Zebrafish",
              initial = TRUE,
              batch_int = FALSE,
              integrate_by = "orig.file", # each channel was separate biological replicate
              harmony = TRUE,
              contamination_threshold = 6, 
              contamination = list(rod = c(NA), cone = c(NA), hc = c(NA), bc = c(NA), ac = NA, rgc = NA, mg = c(0:9, 12, 13, 18, 19, 21, 23, 24)), 
              nFeature_threshold = -1.3,
              doublet_finder = FALSE,
              manual_annotation = list(SAC = c(7,13), A17 = NA, nNOS = c(15), CA1 = c(5), CA2 = c(NA), VG1 = c(NA), A2 = c(4), SEG = c(0), VG3 = NA), 
              sac_annotation = FALSE,
              de_expression = TRUE, 
              save = FALSE)

render_report(species = "Goldfish",
              initial = FALSE,
              batch_int = TRUE,
              integrate_by = "animal", 
              harmony = FALSE,
              contamination_threshold = 6, 
              contamination = list(rod = NA, cone = c(NA), hc = c(NA), bc = NA, ac = NA, rgc = NA, mg = NA), 
              nFeature_threshold = -3,
              doublet_finder = TRUE,
              manual_annotation = list(SAC = c(1), A17 = c(15), nNOS = c(4,5,16), CA1 = c(19), CA2 = c(4), VG1 = c(NA), A2 = c(NA), SEG = c(NA), VG3 = 10), 
              sac_annotation = FALSE,
              de_expression = TRUE, 
              save = TRUE)

render_report(species = "Lamprey",
              initial = TRUE,
              batch_int = TRUE,
              integrate_by = "animal",
              harmony = TRUE,
              contamination_threshold = 6, 
              contamination = list(rod = NA, cone = c(NA), hc = c(NA), bc = NA, ac = NA, rgc = NA, mg = NA), 
              nFeature_threshold = -3,
              doublet_finder = TRUE,
              manual_annotation = list(SAC = c(0,10,32), A17 = NA, nNOS = c(NA), CA1 = c(28), CA2 = c(31), VG1 = NA, A2 = NA, SEG = c(NA), VG3 = 3), 
              sac_annotation = FALSE,
              de_expression = TRUE, 
              save = TRUE)



# Literature proportions
rmarkdown::render("3_literature_proportions.Rmd", output_file = "html_reports/LiteratureProportionAnalysis.html")

# SAC sub-clustering
rmarkdown::render("7_sac_analysis.Rmd", output_file = "html_reports/SacSubclustering_v4.html")

# Cluster reproducibility for each species 
apply(subset(parameters, species %in% c("Opossum", "Goldfish", "Lamprey")), 1, function(row) {
  rmarkdown::render("4_cluster_reproducibility.Rmd", 
                    output_file = paste0("html_reports/ClusterReproducibility-", row["species"], ".html"), 
                    params = list(filepath = paste0("../../Species_Objects/", row["species"], "AC_v5.rds"), 
                                  nPermutations = 100, 
                                  method = row["method"], 
                                  group.by = row["group.by"], 
                                  run.perm = TRUE))
})


# OrthoType analysis for vertebrates
rmarkdown::render("2_orthotype_analysis.Rmd", 
                  output_file = "html_reports/OrthoTypeAC-Vertebrates_v3.html", 
                  params = list(biomart = FALSE,
                                prep_ortho = FALSE,
                                batch_int =  FALSE,
                                harmony =  FALSE,
                                save = FALSE,
                                bigmatrix =  FALSE,
                                shuffling =  FALSE,
                                StopCheckPoint1 =  FALSE,
                                StopCheckPoint2 =  FALSE,
                                StopCheckPoint3 =  FALSE,
                                expr_corr =  FALSE,
                                FindCorrectedGenes = FALSE,
                                merge = FALSE))

# Cluster reproducibility for OT analyses
rmarkdown::render("4_cluster_reproducibility.Rmd", 
                  output_file = "html_reports/ClusterReproducibility-VertebrateAC.html", 
                  params = list(filepath = "../../Ortho_Objects/vertebrateAC_v2_100.rds", 
                                nPermutations = 10, 
                                method = "seurat", 
                                group.by = "species", 
                                run.perm = TRUE))

# Run the TF analysis
rmarkdown::render("5_tf_analysis.Rmd")

# Run 6_samap_analysis.ipynb in jupyter or VSCode

# Run the starburst amacrine cell analysis
rmarkdown::render("7_sac_analysis.Rmd")

# Make all the figures
rmarkdown::render("8_figures.Rmd")
