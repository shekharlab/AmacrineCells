from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                            sankey_plot, chord_plot, CellTypeTriangles, 
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles)
import samap.utils
from samalg import SAM
import pandas as pd
import scanpy as sc
import numpy as np
import seaborn as sn 
import matplotlib.pyplot as plt
from scipy.stats import zscore
from scipy.stats import pearsonr
from sklearn.metrics import mean_squared_error
from sklearn.manifold import MDS 
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from functools import reduce
import os
import glob
import anndata as ad
from multiprocessing import Pool
from itertools import repeat
import pickle

def SaveTables(results, prefix):
    [results[i][1].to_csv(prefix + '_' + str(i) + '.csv') for i in range(len(results))]

def MeanTable(results, **kwargs):
    # element-wise mean
    plt.rcParams['figure.figsize'] = [12, 10]
    tables = [result[1] for result in results]
    mean = pd.DataFrame(np.mean(tables, axis=0), columns=tables[0].columns, index=tables[0].index)
    PlotHeatmap(mean, **kwargs)
    
def StdTable(results):
    # element-wise std
    plt.rcParams['figure.figsize'] = [12, 10]
    tables = [result[1] for result in results]
    mean = pd.DataFrame(np.std(tables, axis=0), columns=tables[0].columns, index=tables[0].index)
    PlotHeatmap(mean)

# the score
def ScoreAlignment(df):

    df_nan = df.copy(deep = True)

    # Make interspecies comparisons NA
    for row in df_nan:
        for column in df_nan:
            if row[:2] == column[:2]:
                df_nan.loc[row,column] = np.nan

    # Rearraging Data
#     new_order = [x for x in new_order if x in df_nan.columns.values]

    # Reorder the rows and columns
#     df_reordered = df_nan.reindex(new_order)[new_order]
    df_reordered = df_nan

    # Metrics calculation
    goodsum_ = []
    badsum_ = []
    for row in df_reordered.index:
        for column in df_reordered.columns:
            if row[3:] == column[3:]:
                value1 = df_reordered.loc[row, column]
                if not pd.isna(value1):  # Check if value is not NaN
                    goodsum_.append(value1)  # Add value to the sum
            else:
                value2 = df_reordered.loc[row, column]
                if not pd.isna(value2):  
                    badsum_.append(value2) 

    min_ = -1*len(badsum_)
    max_ = 1*len(goodsum_)
    # score = sum(goodsum_) - sum(badsum_)
    
    # New metric ranges from -1 to 1 where 0 means similar enrichment of both good and bad edge weights
    accuracy = np.mean(goodsum_) - np.mean(badsum_) # old score: accuracy = (score - min_)/(max_ - min_)
    
    return accuracy, sum(goodsum_)/max_, -sum(badsum_)/min_

def unique(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

def SAMapTrialDebug(filenames, gnnm = None, color = 'cell_class', plot = False, types_remove = None, remove_genes = None, blast_dir = None, keys = None, **kwargs):

    print(filenames)
    identlist = list(filenames.keys())
    
    # Load objects
    print('Loading objects...')
    objects = [ad.read_h5ad(file) for species, file in filenames.items()]
    objects

    # Remove types; doing this in Seurat before export
    if types_remove is not None:
        print('Removing types...')
        objects2 = [SubsetAnndata(object, types_remove = types_remove) for object in objects]
        objects = objects2.copy()
    
    # Remove genes
    if remove_genes is not None:
        print('Removing genes...')
        objects2 = []
        for adata in objects:
            other_genes = [name for name in adata.var_names if not name in remove_genes]
            objects2.append(adata[:, other_genes])
        objects = objects2.copy()
    
    # Take a sample
    print('Downsampling...')
    dsList = [DownsampleAnndata(object, **kwargs) for object in objects]
    
    # Run SAM
    samList = [SamPreprocessing(object) for object in dsList]
    samapDict = dict((identlist[i], samList[i]) for i, species in enumerate(filenames))
    samapDict
    
    # Build SAMap object
    if gnnm is not None:
        print('Using specified gnnm as input...')
        sm = RunSAM(samapDict, gnnm = gnnm)
    else:
        print('Computing gnnm from blast dir...')
        sm = RunSAM(samapDict, blast_dir = blast_dir, keys = keys)
    
    # Run SAMap
    sm, MappingTable, gpf = RunSAMap(sm, 
                                     NUMITERS = 3, 
                                     neigh_from_keys = dict((species, True) for species in identlist))
    return sm, MappingTable, gpf


def SaveSamapData(res, prefix, out_dir = 'samap/', save_gene_pairs = False, align_thr = 0.2, 
                  object_dir = None, species = None, n = None, thr = 0.05):
    """
    Save SAMap analysis results, including alignment tables, clustering metadata, and UMAP embeddings.

    This function exports the results of a SAMap run to disk in multiple formats for downstream analysis.
    It saves the SAMap object itself, the mapping table for R analysis, and optionally gene pair alignments.
    Additionally, it performs Leiden clustering at multiple resolutions and exports UMAP coordinates. Note: you
    will get a division by zero error if you try to use this function and have only one value for "annotated" in 
    any species. You need multiple cell types in order to find gene pairs!

    Parameters
    ----------
    res : tuple
        A tuple containing the SAMap results:
        - `res[0]`: SAMap object.
        - `res[1]`: Mapping table (pandas DataFrame).
        - `res[2]`: Gene pair finder object (optional, used if `save_gene_pairs=True`).
    prefix : str
        Prefix for output filenames (used to distinguish datasets or experiments).
    out_dir : str, optional
        Output directory for CSV exports. Default is `'samap/'`.
    save_gene_pairs : bool, optional
        Whether to compute and save aligned gene pairs. Default is `False`.
    align_thr : float, optional
        Alignment threshold passed to `find_all()` when saving gene pairs. Default is `0.2`.

    Outputs
    -------
    Files are written to disk:
    - `{out_dir}/{prefix}-mappingtable.csv` : Cross-species mapping table.
    - `{out_dir}/{prefix}-umap.csv` : Cell-level UMAP coordinates and Leiden cluster assignments.
    - `"../../../storage/samap/pkl/{prefix}-sm.samap"` : Serialized SAMap object.
    - `{out_dir}/{prefix}-genepairs.csv` : (optional) Gene pairs above the given threshold.

    Notes
    -----
    This function performs Leiden clustering at multiple resolutions (0.5 to 5 in steps of 0.5)
    and stores the results in `samap.adata.obs` before exporting the combined metadata.
    """

    # Export alignment scores to csv for R
    res[1].to_csv(out_dir + prefix + '-mappingtable.csv')

    # Save samap object
    if object_dir is not None: 
        samap.utils.save_samap(res[0], object_dir + '/' + prefix + "-sm.samap")

    # UMAP coordinates and labels
    res[0].samap.leiden_clustering(res=0.5)
    res[0].samap.adata.obs['leiden_clusters_0.5'] = res[0].samap.adata.obs['leiden_clusters']
    res[0].samap.leiden_clustering(res=1)
    res[0].samap.adata.obs['leiden_clusters_1'] = res[0].samap.adata.obs['leiden_clusters']
    res[0].samap.leiden_clustering(res=1.5)
    res[0].samap.adata.obs['leiden_clusters_1.5'] = res[0].samap.adata.obs['leiden_clusters']
    res[0].samap.leiden_clustering(res=2)
    res[0].samap.adata.obs['leiden_clusters_2'] = res[0].samap.adata.obs['leiden_clusters']
    res[0].samap.leiden_clustering(res=3)
    res[0].samap.adata.obs['leiden_clusters_3'] = res[0].samap.adata.obs['leiden_clusters']
    res[0].samap.leiden_clustering(res=3.5)
    res[0].samap.adata.obs['leiden_clusters_3.5'] = res[0].samap.adata.obs['leiden_clusters']
    res[0].samap.leiden_clustering(res=4)
    res[0].samap.adata.obs['leiden_clusters_4'] = res[0].samap.adata.obs['leiden_clusters']
    res[0].samap.leiden_clustering(res=4.5)
    res[0].samap.adata.obs['leiden_clusters_4.5'] = res[0].samap.adata.obs['leiden_clusters']
    res[0].samap.leiden_clustering(res=5)
    res[0].samap.adata.obs['leiden_clusters_5'] = res[0].samap.adata.obs['leiden_clusters']
    df = res[0].samap.adata.obs.copy()
    df["UMAP_1"] = res[0].samap.adata.obsm["X_umap"][:, 0]
    df["UMAP_2"] = res[0].samap.adata.obsm["X_umap"][:, 1]
    df.to_csv(out_dir + prefix + "-umap.csv", index=True)  # index=True keeps cell IDs as first column

    # Gene pairs
    if save_gene_pairs:
        gpf = res[2]
        if species is None: 
            gene_pairs = gpf.find_all(align_thr=align_thr, w1t = 0, w2t = 0, n = n, thr = thr)
        else: 
            gene_pairs = gpf.find_all_species(align_thr=align_thr, species = species, w1t = 0, w2t = 0, n = n, thr = thr)
        gene_pairs.to_csv(out_dir + prefix + '-genepairs.csv')


def SAMapTrial2(filenames, gnnm = None, color = 'cell_class', plot = False, types_remove = None, 
                class_remove = None, pairwise = True, remove_genes = None, blast_dir = None, keys = None, 
                type_downsample = 100, class_downsample = 100000, seed = 12345, verbose = False):

    print(filenames)
    identlist = list(filenames.keys())
    
    # Load objects
    print('Loading objects...')
    objects = [ad.read_h5ad(file) for species, file in filenames.items()]
    objects

    # Remove cell classes; doing this in Seurat before export
    if class_remove is not None:
        print('Removing classes...')
        objects2 = [SubsetAnndata(object, types_remove = class_remove, cluster_key = 'cell_class') for object in objects]
        objects = objects2.copy()

    # Remove types; doing this in Seurat before export
    if types_remove is not None:
        print('Removing types...')
        objects2 = [SubsetAnndata(object, types_remove = types_remove) for object in objects]
        objects = objects2.copy()
    
    # Remove genes
    if remove_genes is not None:
        print('Removing genes...')
        objects2 = []
        for adata in objects:
            other_genes = [name for name in adata.var_names if not name in remove_genes]
            objects2.append(adata[:, other_genes])
        objects = objects2.copy()
    
    # Downsample as type level
    print('Downsampling types...')
    dsList1 = [DownsampleAnndata(object, cluster_key = 'annotated', downsample = type_downsample, seed = seed, verbose = verbose) for object in objects]

    # Downsample at class level
    print('Downsampling class...')
    dsList = [DownsampleAnndata(object, cluster_key = 'cell_class2', downsample = class_downsample, seed = seed, verbose = verbose) for object in dsList1]
    
    # Run SAM
    samList = [SamPreprocessing(object) for object in dsList]
    samapDict = dict((identlist[i], samList[i]) for i, species in enumerate(filenames))
    
    # Build SAMap object
    if gnnm is not None:
        print('Using specified gnnm as input...')
        sm = RunSAM(samapDict, gnnm = gnnm)
    else:
        print('Computing gnnm from blast dir...')
        sm = RunSAM(samapDict, blast_dir = blast_dir, keys = keys)
    
    # Run SAMap
    sm, MappingTable, gpf = RunSAMap(sm, 
                                     NUMITERS = 3, 
                                     pairwise = pairwise,
                                     neigh_from_keys = dict((species, True) for species in identlist))

    
    # Add annotated
    num_species = len(sm.ids)
    annotation_column = "annotated;"*(num_species-1) + "annotated_mapping_scores"
    sm.samap.adata.obs['annotated'] = sm.samap.adata.obs[annotation_column]
    # sm.samap.adata.obs['species_annotated'] = [i.split('_', 1)[1] for i in sm.samap.adata.obs[annotation_column].tolist()]
    sm.samap.adata.obs['cell_class'] = [x.split('-')[0].split('_')[1] for x in sm.samap.adata.obs['annotated']]
    sm.samap.adata.obs['cell_class2'] = np.concatenate([sm.sams[id].adata.obs["cell_class2"] for id in sm.sams.keys()])
    
    if plot:
        plt.rcParams['figure.figsize'] = [5, 5]
        PlotUMAP2(sm.samap.adata, color = ['species', color])

        # plt.rcParams['figure.figsize'] = [5, 5]
        # PlotUMAP2(sm, color = color)
        
    return sm, MappingTable, gpf

def SamPreprocessing(adata):
    sam=SAM(counts = adata)
    sam.preprocess_data(
                        sum_norm="cell_median",
                        norm="log",
                        thresh_low=0.0,
                        thresh_high=0.96, # Key parameter (this was set to 0.96 by default SAMap)
                        min_expression=1,
                    )
    sam.run(
                        preprocessing="StandardScaler",
                        npcs=100,
                        weight_PCs=False,
                        k=20,
                        n_genes=3000,
                        weight_mode='rms',
                        projection="none" # turn off umap calculation
                    )
    return sam

def SubsetAnndata(adata, types_remove, cluster_key = 'annotated'):
    adatas = [adata[adata.obs[cluster_key] == clust] for clust in adata.obs[cluster_key].astype('category').cat.categories if clust not in types_remove]
    adata_merge = adatas[0].concatenate(*adatas[1:])
    return adata_merge

def DownsampleAnndata(adata, downsample = 100, cluster_key = "annotated", seed = 1, verbose = False):

    # Make field categorical 
    adatas = [adata[adata.obs[cluster_key] == clust] for clust in adata.obs[cluster_key].astype('category').cat.categories]

    # Downsample if there's more cells than downsample 
    for dat in adatas:
        if dat.n_obs > downsample:
             sc.pp.subsample(dat, n_obs=downsample, random_state=seed)

    # Merge
    adata_downsampled = adatas[0].concatenate(*adatas[1:])

    # Print numbers
    # for adata in adata_downsampled:
        # print(f"\nSpecies: {adata.obs['species'].iloc[0]}")
    if verbose: 
        print(adata_downsampled.obs[cluster_key].value_counts())
    
    return adata_downsampled

def SamapBootstrap(filenames, gnnm, save = None):
    
    # delete preprocessed files if they exist
    for ident, file in filenames.items(): 
        file_to_delete = file.split('.h5ad')[0]+'_pr.h5ad'
        if os.path.isfile(file_to_delete):
            os.remove(file_to_delete)
            print('Removed ' + file_to_delete)
    
    sm = RunSAM(filenames.copy(), 
                gnnm = gnnm)
    
    sm, MappingTable, gpf = RunSAMap(sm, 
                                     NUMITERS = 3, 
                                     neigh_from_keys = dict((species, True) for species in filenames.keys()))#{'ze':True, 'ch':True, 'li':True, 'op':True, 'rn':True, 'hs':True})

    if save is not None:
        samap.utils.save_samap(sm, save)

    plt.rcParams['figure.figsize'] = [12, 10]
    PlotHeatmap(MappingTable, sm)
#     gene_pairs = gpf.find_all(align_thr=0.10)
#     print(gene_pairs)
    
    return sm, MappingTable, gpf

def RunSAM(filenames, keys = None, blast_dir = 'maps/', gnnm = None):

    if keys is None:
        sm = SAMAP(
                    filenames,
                    f_maps = blast_dir,
                    keys = dict((species, 'annotated') for species in filenames.keys()),
                    gnnm = gnnm, 
                    save_processed=False #if False, do not save the processed results to `*_pr.h5ad`
                )
    else: 
        sm = SAMAP(
                    filenames,
                    f_maps = blast_dir,
                    names = keys,
                    keys = dict((species, 'annotated') for species in filenames.keys()),
                    save_processed=False #if False, do not save the processed results to `*_pr.h5ad`
                )
    return sm

def RunSAM2(filenames, keys = None, blast_dir = 'maps/', gnnm = None):

    sams = []
    for name, file in filenames.items():
        sam=SAM()
        print("Loading ", file)
        sam.load_data(file)
        sam.preprocess_data(
                            sum_norm="cell_median",
                            norm="log",
                            thresh_low=0.0,
                            thresh_high=0.96, # Key parameter (this was set to 0.96 by default SAMap)
                            min_expression=1,
                        )
        sam.run(
                            preprocessing="StandardScaler",
                            npcs=100,
                            weight_PCs=False,
                            k=20,
                            n_genes=3000,
                            weight_mode='rms'
                        )
        sams.append(sam)

    print(dict(zip(filenames.keys(), sams)))
    
    sm = SAMAP(
            dict(zip(filenames.keys(), sams)),
            f_maps = blast_dir,
            gnnm = gnnm, 
#             names = keys,
            keys = dict((species, 'annotated') for species in filenames.keys())
        )

    return sm

def RunSAMap(sm, n_top = 0, **kwargs):
    
    # Run SAMap
    sm.run(**kwargs) #pairwise=pairwise, NUMITERS = NUMITERS)

    # Mapping
    keys = dict((species, 'annotated') for species in sm.ids)
    D,MappingTable = get_mapping_scores(sm, keys, n_top = n_top)
    
    gpf = GenePairFinder(sm,keys=keys)
    
    return sm, MappingTable, gpf

def MapTypes(filenames, keys = None, NUMITERS=3):
    
    # Run SAM
    sm = RunSAM(filenames, keys)
    
    # Run SAMap    
    return RunSAMap(sm, NUMITERS=NUMITERS)

def grouped_obs_mean(adata, group_key, layer=None, gene_symbols=None):
    if layer is not None:
        getX = lambda x: x.layers[layer]
    else:
        getX = lambda x: x.X
    if gene_symbols is not None:
        new_idx = adata.var[idx]
    else:
        new_idx = adata.var_names

    grouped = adata.obs.groupby(group_key)
    out = pd.DataFrame(
        np.zeros((adata.shape[1], len(grouped)), dtype=np.float64),
        columns=list(grouped.groups.keys()),
        index=adata.var_names
    )

    for group, idx in grouped.indices.items():
        X = getX(adata[idx])
        out[group] = np.ravel(X.mean(axis=0, dtype=np.float64))
    return out

def RowNorm(dataframe):
    return dataframe.div(dataframe.sum(axis=1), axis=0)
    
def translate_feature_space(sm, ref, pred, average = True, scale = True, 
                            rownorm = False, scale_prediction = True, gnnm_refined = None):
    
    # Get gene names
    genes1 = [x for x in sm.gns if x.startswith(ref)]
    genes2 = [x for x in sm.gns if x.startswith(pred)]
    
    # Get average expression
    if average:
        X = grouped_obs_mean(sm.samap.adata[sm.samap.adata.obs['species'].isin([ref])][:, genes1], "annotated")
        Y = grouped_obs_mean(sm.samap.adata[sm.samap.adata.obs['species'].isin([pred])][:, genes2], "annotated")
    else: 
        X_adata = sm.samap.adata[sm.samap.adata.obs['species'].isin([ref])][:, genes1]
        X = pd.DataFrame(X_adata.X.toarray().T, index=genes1, columns=X_adata.obs.index)
        print(X.shape)
        Y_adata = sm.samap.adata[sm.samap.adata.obs['species'].isin([pred])][:, genes2]
        Y = pd.DataFrame(Y_adata.X.toarray().T, index=genes2, columns=Y_adata.obs.index)
        print(Y.shape)
    
    # Scale 
    if scale:
        X = pd.DataFrame(scaler.fit_transform(X.T).T, columns=X.columns, index=X.index)
#         X = X.apply(zscore, axis=1)
        X = X.fillna(0)
        Y = pd.DataFrame(scaler.fit_transform(Y.T).T, columns=Y.columns, index=Y.index)
#         Y = Y.apply(zscore, axis=1)
        Y = Y.fillna(0)
    
    # Generate prediction based on reference species and SAMap homology matrix
    if gnnm_refined is None:
        A = sm.gnnm_refined[np.ix_(np.flatnonzero(np.char.startswith(sm.gns, pred)),
                                   np.flatnonzero(np.char.startswith(sm.gns, ref)))]
    else:
        A = gnnm_refined[np.ix_(np.flatnonzero(np.char.startswith(sm.gns, pred)),
                                np.flatnonzero(np.char.startswith(sm.gns, ref)))]
        
    homology_matrix = pd.DataFrame(A.toarray(), 
                                   columns=genes1, 
                                   index=genes2)
    
    # Normalize homology matrix so that rows sum to 1 or cols sum to 1
    if(rownorm): 
        homology_matrix = RowNorm(homology_matrix)
    else: 
        homology_matrix = RowNorm(homology_matrix.T).T
    
    homology_matrix = homology_matrix.fillna(0)
    
    prediction = pd.DataFrame(homology_matrix @ X, 
                              index=genes2)
    
    # Scale prediction
    if(scale_prediction):
        prediction = pd.DataFrame(scaler.fit_transform(prediction.T).T, 
                                  columns=prediction.columns, 
                                  index=prediction.index)
        prediction = prediction.fillna(0)

    # Check that X is in right order
    if X.index.to_list() == genes1:
        print("The lists are identical")
    else:
        print("The lists are not identical")
    
    return X, Y, prediction, homology_matrix

def standardize_expression_matrices(pd_dict, ids, use_intersection = True):
    
    
    for k,v in enumerate(pd_dict):
        
        # Copy
        pd_dict[k] = pd_dict[k].copy(deep=True)
        
        # Add species name to each df
        pd_dict[k].columns = ids[k] + '_' + pd_dict[k].columns
        
        # Remove genes with zero variance
        pd_dict[k] = pd_dict[k].loc[~(pd_dict[k]==0).all(axis=1)]

    # Merge by intersection
    if use_intersection:
        merged = reduce(lambda x, y: pd.merge(x, y, how = 'inner', left_index=True, right_index=True), pd_dict)
    else:
        merged = pd.concat(pd_dict, axis=1)
        merged = merged.fillna(0)
    
    # Print
    print(merged.head())
    print('Shape: ', merged.shape)
    
    return merged

def linear_combination(sm, ref, pred, scale = False, normalize = True, use_A = False, test = False):
    
    # Get gene names
    genes1 = [x for x in sm.gns if x.startswith(ref)]
    genes2 = [x for x in sm.gns if x.startswith(pred)]
    
    # Get average expression
    X = grouped_obs_mean(sm.samap.adata[sm.samap.adata.obs['species'].isin([ref])][:, genes1], "annotated")
    Y = grouped_obs_mean(sm.samap.adata[sm.samap.adata.obs['species'].isin([pred])][:, genes2], "annotated")
    
    if scale:
        X = X.apply(zscore, axis=1)
        X = X.fillna(0)
        Y = Y.apply(zscore, axis=1)
        Y = Y.fillna(0)
    
    # Generate prediction based on reference species and SAMap homology matrix
    A = sm.gnnm_refined[np.ix_(np.flatnonzero(np.char.startswith(sm.gns, pred)),
                               np.flatnonzero(np.char.startswith(sm.gns, ref)))]
    homology_matrix = pd.DataFrame(A.toarray(), columns=genes1, index=genes2)
    
    # Check that X is in right order
    if X.index.to_list() == genes1:
        print("The lists are identical")
    else:
        print("The lists are not identical")
    
    if use_A:
        prediction = pd.DataFrame(A @ X, index=genes2, columns=X.columns)
    else: 
        prediction = pd.DataFrame(homology_matrix @ X, index=genes2)

    if test: 
        return prediction
    
    prediction['bias'] = 1
    
    # Apply linear regression
    w = np.linalg.inv(np.transpose(prediction) @ (prediction)) @ (np.transpose(prediction) @ Y)
    return prediction, Y
    #row_names = dict((idx, ele) for idx,ele in enumerate(w.columns))
    row_names = dict((idx, ele) for idx,ele in enumerate(X.columns.values))
    w = w.rename(index = row_names)
    
    # Normalization (each set of weights will sum to 1)
    if normalize:
#         w[w.columns] = w[w.columns] / w[w.columns].sum()
        wnorm=(w-w.min())/(w.max()-w.min()) # Min max normalization
    
    # Heatmap
    ax = sn.heatmap(data=wnorm, 
                   annot=True, 
                   annot_kws={"fontsize":8}, 
                   cmap='coolwarm',
                   center=0)

    plt.xlabel(pred) # x-axis label with fontsize 15
    plt.ylabel(ref) # y-axis label with fontsize 15
    plt.show()
    
    return X, Y, w, prediction, homology_matrix

def GoodnessOfFit(w, prediction, Y, pred = "Inferred", ref = "True", rsquared = True):
    best = prediction @ w.to_numpy() 

    # Compute R^2
    if rsquared: 
        R_squared = [[pearsonr(best.iloc[:,i], Y.iloc[:,j]).statistic**2 for i in range(len(best.columns))] for j in range(len(Y.columns))]
    else: 
        R_squared = [[mean_squared_error(Y.iloc[:,j], best.iloc[:,i], squared=False) for i in range(len(best.columns))] for j in range(len(Y.columns))]
    R_squared_matrix = pd.DataFrame(np.array(R_squared))
    
    # Annotate rows and columns
    R_squared_matrix.columns = w.columns
    R_squared_matrix.index = w.columns
    
    # Heatmap
    ax = sn.heatmap(data=R_squared_matrix, 
               annot=True, 
               annot_kws={"fontsize":8},
               cmap='coolwarm',
              center = np.min(R_squared_matrix))
    
    plt.xlabel(pred) # x-axis label with fontsize 15
    plt.ylabel(ref) # y-axis label with fontsize 15
    
    return R_squared_matrix

def GeneSimHeatmap(gpf, row_normalize = False):

    gene_pairs = gpf.find_all(align_thr=0)
    list_similarities = gene_pairs.loc[:, ~gene_pairs.columns.str.endswith(('_pval1', '_pval2'))].count()

    gene_summary = pd.DataFrame({'type1': [i.split(';', 1)[0] for i in list_similarities.index], 
                                  'type2': [i.split(';', 1)[1] for i in list_similarities.index],
                                  'n_genes':list_similarities})

    # make unique, sorted, common index
    idx = sorted(set(gene_summary['type1']).union(gene_summary['type2']))

    # reshape
    gene_matrix = (gene_summary.pivot(index='type1', columns='type2', values='n_genes')
       .reindex(index=idx, columns=idx)
       .fillna(0, downcast='infer')
       .pipe(lambda x: x+x.values.T)
     )

    # Row normalize
    if row_normalize:
        row_sums = gene_matrix.sum(axis=1)
        gene_matrix = gene_matrix / row_sums[:, np.newaxis]

    sn.heatmap(data=gene_matrix, 
               annot=True, 
               annot_kws={"fontsize":6}) 
    
    return gene_matrix


def PlotPCA(data, color = True, binary = None):
    # From https://builtin.com/machine-learning/pca-in-python
    # For subplots https://stackoverflow.com/questions/20073017/return-a-subplot-from-a-function

    X = data.copy(deep = True)
    if binary is not None:
        X[data >= binary] = 1
        X[data < binary] = 0
    
    # Standardizing the features
    x = StandardScaler().fit_transform(X.T)

    # PCA
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data = principalComponents
                 , columns = ['principal component 1', 'principal component 2'])
    principalDf.shape
    celltypes = X.columns.values.tolist()
    principalDf['celltype'] = celltypes
    finalDf = principalDf
    print(finalDf.head())
    # finalDf = pd.concat([principalDf, pd.DataFrame(target = X.columns.values.tolist())], axis = 1)
    
    if color:
        colors = [string.split('_')[1] for string in celltypes]
        for index, item in enumerate(colors):
            if item == "UV":
                colors[index] = "purple"
            if item == "rod":
                colors[index] = "grey"
            if item == "accessory":
                colors[index] = "cyan"
            if item == "principle":
                colors[index] = "magenta"
            if item == "principal":
                colors[index] = "magenta"
    else: 
        colors = "black"
            
    # Plot
    fig = plt.figure(figsize = (5,5))
    ax = fig.add_subplot(1,1,1) 
    ax.set_xlabel('PC1 (' + str(round(pca.explained_variance_ratio_[0]*100, 1)) + '%)', fontsize = 10)
    ax.set_ylabel('PC2 (' + str(round(pca.explained_variance_ratio_[1]*100, 1)) + '%)', fontsize = 10)
    ax.set_title('PCA', fontsize = 20)
#     ax.grid()
    ax.scatter(finalDf.loc[:, 'principal component 1'], 
               finalDf.loc[:, 'principal component 2'], 
               c = colors, 
               s = 50)
#     ax.legend(celltypes)

    for i, txt in enumerate(celltypes):
        ax.annotate(txt, (finalDf.loc[i, 'principal component 1'], finalDf.loc[i, 'principal component 2']))
        
    return finalDf
        

def PlotUMAP2(adata, color = 'species', palette = None, size=10, shuffle = True):
    # adata = sm.samap.adata.copy()
    obs = adata.obs.replace('unassigned',np.NaN)

    if shuffle: 
        adata = adata[np.random.permutation(adata.n_obs), :]
    
    sc.pl.umap(adata, 
           size=size,
           color=color,
           palette=palette)
    
def PlotUMAP(sm, size=10, palette = None):
    num_species = len(sm.ids)
    obs = sm.samap.adata.obs.replace('unassigned',np.NaN)
    annotation_column = "annotated;"*(num_species-1) + "annotated_mapping_scores"
    sm.samap.adata.obs['annotated'] = [i.split('_', 1)[1] for i in sm.samap.adata.obs[annotation_column].tolist()]
    sc.pl.umap(sm.samap.adata, 
           size=size,
           color=['species','annotated'],
           palette={'ze':'tab:green', 'op':'tab:blue','ch':'tab:red','li':'yellow','rn':'magenta', 'sq':'magenta','hs':'cyan',
                    'UV':'tab:purple','green':'tab:green','blue':'tab:blue','red':'tab:red','rod':'tab:grey', 
                    'OFF_bipolar':'black', 'ON_bipolar':'brown', 'AC':'black', 'RGC':'brown','HC':'brown', 
                    'SAC':'brown','VG3':'black',
                    'principle':'magenta', 'principal':'magenta', 'accessory':'cyan', 'opn1mw4/opn1lw1':'brown'})
    
def PlotHeatmap(MappingTable, **kwargs):
#     keys = dict((species, 'annotated') for species in sm.ids)
#     D,MappingTable = get_mapping_scores(sm, keys, n_top = 0)

    # Make interspecies comparisons NA
    for row in MappingTable:
        for column in MappingTable:
            if row[:2] == column[:2]:
                MappingTable.loc[row,column] = np.nan

    # Reorder
    ids = unique([column.split('_')[0] for column in MappingTable.columns])
    new_order = [species + '_' + type for type in ['UV', 'blue', 'green', 'opn1mw4/opn1lw1', 'red', 'principle','principal','accessory', 'rod'] for species in ids]
    new_order = [x for x in new_order if x in MappingTable.columns.values]
#     np.fill_diagonal(MappingTable.values, np.NaN)
    MappingTable = MappingTable.reindex(new_order)[new_order].round(2) # Round values
    sn.heatmap(data=MappingTable, 
               annot=True, 
               annot_kws={"fontsize":6}, 
               **kwargs)
    
def PlotMDS(MappingTable):
    mds = MDS(n_components=2, random_state=0) 
  
    # Fit the data to the MDS 
    # object and transform the data 
    X_transformed = pd.DataFrame(mds.fit_transform(MappingTable), columns = ['Comp1', 'Comp2'])
    X_transformed.index = MappingTable.index
    print(X_transformed)
    
    # Plot
    fig = plt.figure(figsize = (8,8))
    ax = fig.add_subplot(1,1,1)
    ax.grid()
    ax.set_xlabel('Component 1', fontsize = 15)
    ax.set_ylabel('Component 2', fontsize = 15)
    ax.set_title('MDS', fontsize = 20)
    
    colors = [string.split('_')[1] for string in MappingTable.index]
    for index, item in enumerate(colors):
        if item == "UV":
            colors[index] = "purple"
        if item == "rod":
            colors[index] = "grey"
        if item == "accessory":
            colors[index] = "cyan"
        if item == "principle":
            colors[index] = "magenta"
    
    ax.scatter(X_transformed.iloc[:,0], 
               X_transformed.iloc[:,1], 
               c = colors, 
               s = 50)
    
    for i, txt in enumerate(X_transformed.index):
        ax.annotate(txt, (X_transformed.loc[txt, 'Comp1'], X_transformed.loc[txt, 'Comp2']))

def PlotConnectivity(sm):
    c_matrix = pd.DataFrame(sm.samap.adata.obsp['connectivities'].A)
    cell_names = sm.samap.adata.obs.index

    # Annotations
    annotations = sm.samap.adata.obs['species'].astype(str) + '_' + sm.samap.adata.obs['annotated'].astype(str)
    new_order = [species + '_' + type for type in ['UV', 'blue', 'green', 'red', 'rod','principle','accessory','SAC','VG3','HC','AC'] for species in ['ze', 'ch', 'li', 'op', 'rn', 'hs']]
    new_order = [x for x in new_order if x in annotations.to_list()]
    annotations = pd.Categorical(annotations, categories=new_order)

    # Change the column names
    c_matrix.columns = annotations

    # Change the row indexes
    c_matrix.index = annotations

    # plt.rcParams['figure.figsize'] = [100, 100]
    plt.figure(figsize=(15, 15), dpi = 300) 

    # sn.heatmap(data=c_matrix.drop(columns=['annotation']), annot=False)
    ax = sn.heatmap(data=c_matrix.iloc[annotations.argsort(),annotations.argsort()], annot=False)
    line_indices = [annotations[annotations.argsort()].tolist().index(type) for type in new_order]
    ax.hlines(line_indices, *ax.get_xlim(), colors="white", linewidth=0.1)
    ax.vlines(line_indices, *ax.get_xlim(), colors="white", linewidth=0.1)
    ax


def find_all_species(self, n=None, align_thr=0.1, n_top=0, species = 'ch', **kwargs):
    """Modified find_all: subsets M.values to rows whose index starts with 'species' ident."""
    _, M = samap.analysis.get_mapping_scores(self.sm, self.keys, n_top=n_top)

    # print(M.iloc[0:10,0:4])
    # Subset rows and columns based on index name prefix
    mask = M.index.str.startswith(tuple(species))
    # print(mask)
    M = M.loc[mask, :]
    # print(M.iloc[0:10,0:4])
    # print(M.shape)

    # ax = samap.analysis.q(M.columns)
    # print(ax)
    ax_x = samap.analysis.q(M.index)
    ax_y = samap.analysis.q(M.columns)

    # Continue as before
    data = M.values.copy()
    data[data < align_thr] = 0
    x, y = data.nonzero()
    ct1, ct2 = ax_x[x], ax_y[y]
    # print(ct1)

    if n is not None:
        f1 = ct1 == n
        f2 = ct2 == n
        f = np.logical_or(f1, f2)
    else:
        f = np.array([True] * ct2.size)

    ct1 = ct1[f]
    ct2 = ct2[f]
    ct1, ct2 = np.unique(np.sort(np.vstack((ct1, ct2)).T, axis=1), axis=0).T
    # print(ct1)

    res = {}
    for i in range(ct1.size):
        a = '_'.join(ct1[i].split('_')[1:])
        b = '_'.join(ct2[i].split('_')[1:])
        print(
            f"Calculating gene pairs for mapping: {ct1[i].split('_')[0]},{a} "
            f"to {ct2[i].split('_')[0]},{b}"
        )
        res[f"{ct1[i]};{ct2[i]}"] = self.find_genes_fixed(ct1[i], ct2[i], **kwargs)

    cols = []
    col_names = []
    for k in res:
        col_names.append(k)
        col_names.append(k + "_pval1")
        col_names.append(k + "_pval2")
        cols.append(res[k][0])
        cols.append(res[k][-2])
        cols.append(res[k][-1])

    res = pd.DataFrame(cols, index=col_names).fillna(np.nan).T
    return res

def _find_link_genes_avg_fixed(self, c1, c2, id1, id2, w1t=0.0, w2t=0.0, expr_thr=0.05):
        mus = self.mus
        stds = self.stds
        sams=self.sams

        keys=self.keys
        sam3=self.s3
        gnnm = self.gnnm
        gns = self.gns
        
        xs = []
        for sid in [id1,id2]:
            xs.append(sams[sid].get_labels(keys[sid]).astype('str').astype('object'))
        x1,x2 = xs
        g1, g2 = gns[np.vstack(gnnm.nonzero())]
        gs1,gs2 = samap.analysis.q([x.split('_')[0] for x in g1]), samap.analysis.q([x.split('_')[0] for x in g2])
        filt = np.logical_and(gs1==id1,gs2==id2)
        g1=g1[filt]
        g2=g2[filt]
        sam1,sam2 = sams[id1],sams[id2]
        mu1,std1,mu2,std2 = mus[id1][g1].values,stds[id1][g1].values,mus[id2][g2].values,stds[id2][g2].values

        X1 = samap.analysis._sparse_sub_standardize(sam1.adata[:, g1].X[x1 == c1, :], mu1, std1)
        X2 = samap.analysis._sparse_sub_standardize(sam2.adata[:, g2].X[x2 == c2, :], mu2, std2)
        a, b = sam3.adata.obsp["connectivities"][sam3.adata.obs['species']==id1,:][:,sam3.adata.obs['species']==id2][
            x1 == c1, :][:, x2 == c2].nonzero()
        c, d = sam3.adata.obsp["connectivities"][sam3.adata.obs['species']==id2,:][:,sam3.adata.obs['species']==id1][
            x2 == c2, :][:, x1 == c1].nonzero()            

        pairs = np.unique(np.vstack((np.vstack((a, b)).T, np.vstack((d, c)).T)), axis=0)

        av1 = X1[np.unique(pairs[:, 0]), :].mean(0).A.flatten()
        av2 = X2[np.unique(pairs[:, 1]), :].mean(0).A.flatten()
        sav1 = (av1 - av1.mean()) / av1.std()
        sav2 = (av2 - av2.mean()) / av2.std()
        sav1[sav1 < 0] = 0
        sav2[sav2 < 0] = 0
        val = sav1 * sav2 / sav1.size
        X1.data[:] = 1
        X2.data[:] = 1
        min_expr = (X1.mean(0).A.flatten() > expr_thr) * (
            X2.mean(0).A.flatten() > expr_thr
        )

        w1 = sam1.adata.var["weights"][g1].values.copy()
        w2 = sam2.adata.var["weights"][g2].values.copy()
        w1[w1 < w1t] = 0
        w2[w2 < w2t] = 0
        w1[w1 > 0] = 1
        w2[w2 > 0] = 1
        return val * w1 * w2 * min_expr, samap.utils.to_vn(np.array([g1,g2]).T)


def find_genes_fixed(
        self,
        n1,
        n2,
        w1t=0.2,
        w2t=0.2,
        n_genes=1000,
        thr=1e-2,
    ):
        """Find enriched gene pairs in a particular pair of cell types.
        
        n1: str, cell type ID from species 1
        
        n2: str, cell type ID from species 2
        
        w1t & w2t: float, optional, default 0.2
            SAM weight threshold for species 1 and 2. Genes with below this threshold will not be
            included in any enriched gene pairs.
        
        n_genes: int, optional, default 1000
            Takes the top 1000 ranked gene pairs before filtering based on differential expressivity and
            SAM weights.
        
        thr: float, optional, default 0.01
            Excludes genes with greater than 0.01 differential expression p-value.
            
        Returns
        -------
        G - Enriched gene pairs
        G1 - Genes from species 1 involved in enriched gene pairs
        G2 - Genes from species 2 involved in enriched gene pairs
        pvals1 - pvalues for genes from species 1 involved in enriched gene pairs
        pvals2 - pvalues for genes from species 2 involved in enriched gene pairs
        """
        n1 = str(n1)
        n2 = str(n2)
        id1,id2 = n1.split('_')[0],n2.split('_')[0]
        sam1,sam2=self.sams[id1],self.sams[id2]

        n1,n2 = '_'.join(n1.split('_')[1:]),'_'.join(n2.split('_')[1:])
        assert n1 in samap.analysis.q(self.sams[id1].adata.obs[self.keys[id1]])
        assert n2 in samap.analysis.q(self.sams[id2].adata.obs[self.keys[id2]])
        
        m,gpairs = self._find_link_genes_avg_fixed(n1, n2, id1,id2, w1t=w1t, w2t=w2t, expr_thr=0.05)

        self.gene_pair_scores = pd.Series(index=gpairs, data=m)

        G = samap.analysis.q(gpairs[np.argsort(-m)[:n_genes]])
        G1 = samap.utils.substr(G, ";", 0)
        G2 = samap.utils.substr(G, ";", 1)
        pvals1 = samap.analysis.q(sam1.adata.varm[self.keys[id1] + "_pvals"][n1][G1])
        pvals2 = samap.analysis.q(sam2.adata.varm[self.keys[id2] + "_pvals"][n2][G2])
        filt = np.logical_and(pvals1 < thr,pvals2 < thr)
        G = samap.analysis.q(
            G[filt]
        )
        G1 = samap.utils.substr(G, ";", 0)
        G2 = samap.utils.substr(G, ";", 1)
        _, ix1 = np.unique(G1, return_index=True)
        _, ix2 = np.unique(G2, return_index=True)
        G1 = G1[np.sort(ix1)]
        G2 = G2[np.sort(ix2)]
        return G, G1, G2, pvals1[filt], pvals2[filt]


# Patch it in
GenePairFinder.find_all_species = find_all_species
GenePairFinder.find_genes_fixed = find_genes_fixed
GenePairFinder._find_link_genes_avg_fixed = _find_link_genes_avg_fixed

