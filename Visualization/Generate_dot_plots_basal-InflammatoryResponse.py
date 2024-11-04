#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import os
import ast
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# In[28]:


# # Load the metadata and counts data
#Base directory
base_dir = "/Users/weiwu2/Library/CloudStorage/Box-Box/Research/SingleCell_network/Akshat/run_Pyscenic"
manuscript_dir = "/Users/weiwu2/Library/CloudStorage/Box-Box/Research/SingleCell_network/Manuscript/"
DEG_regulon_dir = f"{base_dir}/results/Cluster_Analysis/Results_061424/DEG_results/DEG_Regulons"

#code directory
code_dir = f"{base_dir}/code/Wei/Run_dot_plot"

#Enrichment_results_dir
enrichment_results_path = f'{code_dir}/Enrichment_results.xlsx'

# # Define the directory where you want to save the PDF
#save_dir = f"{base_dir}/results/Cluster_Analysis/Results_061424/Dot_plots/Basal/Latest/"
save_dir = f"{manuscript_dir}/Figures_Tables/Dot_plots/Latest_v3/"
stats_dir = f"{base_dir}/results/Cluster_Analysis/Results_061424/DEG_results/DEG_Regulons/Basal/no_pval_cutoff/DEG_TargetGenes_Sig_Regulons_fc15/Stats/"

# metadata and counts data directory
metadata_path = f'{code_dir}/metadata_for_dp_scaled.csv'
counts_path_scaled = f'{code_dir}/counts_data_for_dp_scaled.csv'
#counts_path_unscaled = f'{code_dir}/counts_data_for_dp_scaled.csv'

# Load the metadata and counts data
metadata = pd.read_csv(metadata_path, index_col=0)
counts_scaled = pd.read_csv(counts_path_scaled, index_col=0)
#counts_unscaled = pd.read_csv(counts_path_unscaled, index_col=0)

# Create AnnData object
adata_scaled = ad.AnnData(X=counts_scaled.T, obs=metadata)
#adata_unscaled = ad.AnnData(X=counts_unscaled.T, obs=metadata)



# In[4]:


# Filter for 5ht -scaled
adata_5ht_scaled = adata_scaled[adata_scaled.obs['orig.ident'] == '5Ht']

# Filter for 6ho -scaled
adata_6ho_scaled = adata_scaled[adata_scaled.obs['orig.ident'] == '6Ho']

# Filter for basal -scaled
adata_scaled_basal = adata_scaled[adata_scaled.obs['cluster1'] == 'Mammary epithelial cells-Basal']




# In[7]:


u1_5ht_stats = pd.read_csv(f'{stats_dir}/5ht/U1_wt_5Ht.csv')
u1_5ht_stats                               





# # U1_wt

# In[9]:


#sorts = ['Fraction','Mean']

curr_enrichment_df = pd.read_excel(f'{code_dir}/Enrichment_tables/Inflammatory_response_reg_u1_wt.xlsx',index_col=0)

#for sort in sorts:
pdf_filename = os.path.join(save_dir, f"Inflammatory_response_reg_U1_wt_dp_sort_stats.pdf")
with PdfPages(pdf_filename) as pdf:
    curr_metadata = metadata[metadata['cluster1'] == 'Mammary epithelial cells-Basal']
    curr_sub_info = [clust.split('l_')[1] for clust in curr_metadata['sub_info']]
    curr_metadata['curr_sub_info'] = curr_sub_info
    adata_scaled_basal.obs = curr_metadata

    # Filter for 5ht and 6ho
    adata_5ht_basal = adata_scaled_basal[adata_scaled_basal.obs['orig.ident'] == '5Ht']
    adata_6ho_basal = adata_scaled_basal[adata_scaled_basal.obs['orig.ident'] == '6Ho']

    adata_s1 = adata_5ht_basal[adata_5ht_basal.obs['curr_sub_info']=='S1']
    adata_s2 = adata_5ht_basal[adata_5ht_basal.obs['curr_sub_info']=='S2']
    adata_s3 = adata_5ht_basal[adata_5ht_basal.obs['curr_sub_info']=='S3']
    adata_s4 = adata_5ht_basal[adata_5ht_basal.obs['curr_sub_info']=='S4']
    adata_u1 = adata_5ht_basal[adata_5ht_basal.obs['curr_sub_info']=='U1_wt']
    adatas = [adata_s1,adata_s2,adata_s3,adata_s4,adata_u1]


    for i in range(len(curr_enrichment_df)):
        tf = curr_enrichment_df['Transcription Factor'][i]
        regulon = tf + '(+)'
        genes = curr_enrichment_df['DETG'][i].split('/')
        genes = [gene for gene in genes if gene != '']
        # Ensure the directory exists
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        # Sort by order from stats file

        # Extract the specific cell (which is a string representation of a list)
        target_genes_string = u1_5ht_stats[u1_5ht_stats['Significant_Regulons'] == regulon]['Significant Target Genes (SameCell)'].iloc[0]

        # Convert the string to a list using ast.literal_eval
        target_genes_list = ast.literal_eval(target_genes_string)

        # Step 1: Create a mapping from the first list to its index position
        index_map = {value: index for index, value in enumerate(target_genes_list)}

        # Step 2: Sort the second list based on the index in the first list
        sorted_genes = sorted(genes, key=lambda x: index_map.get(x, float('inf')))


        sorted_genes = [tf] + sorted_genes   

        # Create a PdfPages object
        fig, axs = plt.subplots(1, 2, figsize=(12, 12))
        sc.set_figure_params(scanpy=True, fontsize=15)

        # Comparison across all 5ht clusters - scaled- mean only expressed - FALSE
        sc.pl.dotplot(
            adata_5ht_basal, 
            var_names=sorted_genes,  # Use sorted genes
            groupby='curr_sub_info', 
            standard_scale='var', 
            dot_max=0.8, 
            dot_min=0.1, 
            color_map='Reds', 
            size_title='Fraction of cells %',
            mean_only_expressed=False,
            swap_axes=True,
            ax=axs[0],
            show=False,
            title='5ht'
        )

        # Comparison across all 6ho clusters - scaled- mean only expressed - FALSE
        sc.pl.dotplot(
            adata_6ho_basal, 
            var_names=sorted_genes,  # Use sorted genes
            groupby='curr_sub_info', 
            standard_scale='var', 
            dot_max=0.8, 
            dot_min=0.1, 
            color_map='Reds', 
            size_title='Fraction of cells %',
            mean_only_expressed=False,
            swap_axes=True,
            ax=axs[1],   
            show=False,
            title='6ho'
        )

        # Adjust layout
        plt.tight_layout()

        # Save the figure to the PDF file
        pdf.savefig(fig)  # Save the current figure
        plt.close(fig)  # Close the figure to free up memory
    


# # DEGs 

# In[20]:


u1_imm_resp_diff = [
    "Fcer1g", "Ctss", "Alox5ap", "Lpl", "Ccr2", "Lyn", "Csf1r", "Tnfaip3", "Ctsc", "Gpx1", "Ncf1", "Fcgr3", "Fcgr2b", 
    "Ccr5", "Adam8", "Metrnl", "Grn", "Ninj1", "Nfkbiz", "Tlr2", "Il1b", "Tnfrsf1b", "Nlrp3", "Stap1", "Adora2b", 
    "Plcg2", "Abhd12", "Gpsm3", "Nupr1", "Sirpa", "Nfkbia", "Nfkb1", "Ddt", "Ccl3", "Daglb", "Tnf", "Cd44", "Zbp1", 
    "Zfp36", "Igf1", "Casp1", "Cd47", "Pycard", "Jak2"
]


# # U1_wt

# In[7]:


sorts = ['Fraction','Mean']

for sort in sorts:
    pdf_filename = os.path.join(save_dir, f"Inflammatory_DEG_U1_wt_dp_sort_{sort}.pdf")
    with PdfPages(pdf_filename) as pdf:
        curr_metadata = metadata[metadata['cluster1'] == 'Mammary epithelial cells-Basal']
        curr_sub_info = [clust.split('l_')[1] for clust in curr_metadata['sub_info']]
        curr_metadata['curr_sub_info'] = curr_sub_info
        adata_scaled_basal.obs = curr_metadata

        # Filter for 5ht and 6ho
        adata_5ht_basal = adata_scaled_basal[adata_scaled_basal.obs['orig.ident'] == '5Ht']
        adata_6ho_basal = adata_scaled_basal[adata_scaled_basal.obs['orig.ident'] == '6Ho']

        adata_s1 = adata_5ht_basal[adata_5ht_basal.obs['curr_sub_info']=='S1']
        adata_s2 = adata_5ht_basal[adata_5ht_basal.obs['curr_sub_info']=='S2']
        adata_s3 = adata_5ht_basal[adata_5ht_basal.obs['curr_sub_info']=='S3']
        adata_s4 = adata_5ht_basal[adata_5ht_basal.obs['curr_sub_info']=='S4']
        adata_u1 = adata_5ht_basal[adata_5ht_basal.obs['curr_sub_info']=='U1_wt']
        adatas = [adata_s1,adata_s2,adata_s3,adata_s4,adata_u1]


        #genes = u1_imm_resp_diff
        # Ensure the directory exists
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        if sort == 'Fraction':

            fraction_expressed = (adata_u1[:, u1_imm_resp_diff].X > 0).mean(axis=0)

            # Convert to a pandas Series for sorting
            fraction_expressed_series = pd.Series(fraction_expressed, index=u1_imm_resp_diff)

            # Step 2: Sort genes by fraction of cells
            sorted_genes = fraction_expressed_series.sort_values(ascending=False).index.tolist()

        else:
            gene_exp = []
            for g in u1_imm_resp_diff:
                mean_exp_gene = []
                for data in adatas:
                    mean_exp_gene.append(data[:, [g]].to_df().mean(axis=0)[0])
                mean_exp_gene = (mean_exp_gene - np.min(mean_exp_gene)) / (np.max(mean_exp_gene) - np.min(mean_exp_gene))
                gene_exp.append(mean_exp_gene[-1])#the -1 index represents cluster u1_wt

            # Step 1: Pair the elements together
            paired_lists = list(zip(u1_imm_resp_diff, gene_exp))

            # Step 2: Sort based on the second list (list2), in descending order
            sorted_pairs = sorted(paired_lists, key=lambda x: x[1], reverse=True)

            # Step 3: Unzip the sorted pairs back into two lists
            sorted_list1, sorted_list2 = zip(*sorted_pairs)

            # Convert them back to lists (since zip returns tuples)
            sorted_genes = list(sorted_list1)


        # Create a PdfPages object
        fig, axs = plt.subplots(1, 2, figsize=(10, 10))

        # Comparison across all 5ht clusters - scaled- mean only expressed - FALSE
        sc.pl.dotplot(
            adata_5ht_basal, 
            var_names=sorted_genes,  # Use sorted genes
            groupby='curr_sub_info', 
            standard_scale='var', 
            dot_max=0.8, 
            dot_min=0.1, 
            color_map='Reds', 
            size_title='Fraction of cells %',
            mean_only_expressed=False,
            swap_axes=True,
            ax=axs[0],
            show=False,
            title='5ht'
        )

        # Comparison across all 6ho clusters - scaled- mean only expressed - FALSE
        sc.pl.dotplot(
            adata_6ho_basal, 
            var_names=sorted_genes,  # Use sorted genes
            groupby='curr_sub_info', 
            standard_scale='var', 
            dot_max=0.8, 
            dot_min=0.1, 
            color_map='Reds', 
            size_title='Fraction of cells %',
            mean_only_expressed=False,
            swap_axes=True,
            ax=axs[1],   
            show=False,
            title='6ho'
        )

        # Adjust layout
        plt.tight_layout()

        # Save the figure to the PDF file
        pdf.savefig(fig)  # Save the current figure
        plt.close(fig)  # Close the figure to free up memory

