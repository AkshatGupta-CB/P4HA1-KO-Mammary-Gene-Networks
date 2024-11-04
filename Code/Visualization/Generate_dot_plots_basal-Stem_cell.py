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
save_dir = f"{manuscript_dir}/Figures_Tables/Dot_plots/Latest_v4/"
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


# In[3]:


adata_scaled


# In[4]:


# Filter for 5ht -scaled
adata_5ht_scaled = adata_scaled[adata_scaled.obs['orig.ident'] == '5Ht']

# Filter for 6ho -scaled
adata_6ho_scaled = adata_scaled[adata_scaled.obs['orig.ident'] == '6Ho']

# Filter for basal -scaled
adata_scaled_basal = adata_scaled[adata_scaled.obs['cluster1'] == 'Mammary epithelial cells-Basal']

# # Filter for 5ht -unscaled
# adata_5ht_unscaled = adata_unscaled[adata_unscaled.obs['orig.ident'] == '5Ht']

# # Filter for 6ho -unscaled
# adata_6ho_unscaled = adata_unscaled[adata_unscaled.obs['orig.ident'] == '6Ho']


# In[21]:


# enrichment_results = pd.ExcelFile(enrichment_results_path)
# sheets = enrichment_results.sheet_names[5:]
# sheets


# In[4]:


u1_5ht_stats = pd.read_csv(f'{stats_dir}/5ht/U1_wt_5Ht.csv')
u1_5ht_stats                           





# # U1_wt

# In[36]:


#sorts = ['Fraction','Mean']

curr_enrichment_df = pd.read_csv(f'{code_dir}/Enrichment_tables/Stem_Cell_U1_wt.tsv',sep='\t',index_col=0)
TFs_of_int = ['Ikzf1(+)', 'Fli1(+)', 'Rel(+)', 'Mafb(+)', 'Jund(+)']
curr_enrichment_df = curr_enrichment_df[curr_enrichment_df['Transcription Factor'].isin(TFs_of_int)].reset_index(drop=True)
#for sort in sorts:
pdf_filename = os.path.join(save_dir, f"Stem_Cell_reg_U1_wt_dp_sort_stats_v4.pdf")
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
        regulon = curr_enrichment_df['Transcription Factor'][i]
        tf = regulon.split('(')[0]
        genes = curr_enrichment_df['DETG'][i].split(', ')
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
        if len(sorted_genes) <= 8:
            size = 6
        elif len(sorted_genes) > 8 and len(sorted_genes) <= 10:
            size = 8
        else:
            size = 10

        fig, axs = plt.subplots(1, 2, figsize=(10,size))
        sc.set_figure_params(scanpy=True, fontsize=20)
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
    


# # 5ht_s2

# In[27]:


#Read the stats dataframe
s2_5ht_stats = pd.read_csv(f'{stats_dir}5ht/S2_5Ht.csv')
s2_5ht_stats   


# In[25]:


#sorts = ['Fraction','Mean']
curr_enrichment_df = pd.read_csv(f'{code_dir}/Enrichment_tables/Stem_Cell_5ht_S2.tsv',sep='\t',index_col=0)

#for sort in sorts:
pdf_filename = os.path.join(save_dir, f"Stem_Cell_reg_S2_5ht_dp_sort_stats_v4.pdf")
with PdfPages(pdf_filename) as pdf:
    curr_metadata = metadata[metadata['cluster1'] == 'Mammary epithelial cells-Basal']
    curr_sub_info = [clust.split('l_')[1] for clust in curr_metadata['sub_info']]
    curr_metadata['curr_sub_info'] = curr_sub_info
    adata_scaled_basal.obs = curr_metadata

    # Filter for 5ht and 6ho
    adata_5ht_basal = adata_scaled_basal[adata_scaled_basal.obs['orig.ident'] == '5Ht']
    adata_6ho_basal = adata_scaled_basal[adata_scaled_basal.obs['orig.ident'] == '6Ho']

    # Step 1: Calculate fraction of cells expressing each gene
    # Calculate fraction of cells for each gene in the current set
    S2_basal_adata = adata_5ht_basal[adata_5ht_basal.obs['curr_sub_info']=='S2']

    for i in range(len(curr_enrichment_df)):
        regulon = curr_enrichment_df['Transcription Factor'][i]
        tf = regulon.split('(')[0]
        genes = curr_enrichment_df['DETG'][i].split(', ')
        genes = [gene for gene in genes if gene != '']
        # Ensure the directory exists
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)



        # Sort by order from stats file

        # Extract the specific cell (which is a string representation of a list)
        target_genes_string = s2_5ht_stats[s2_5ht_stats['Significant_Regulons'] == regulon]['Significant Target Genes (5Ht_6Ho)'].iloc[0]

        # Convert the string to a list using ast.literal_eval
        target_genes_list = ast.literal_eval(target_genes_string)

        # Step 1: Create a mapping from the first list to its index position
        index_map = {value: index for index, value in enumerate(target_genes_list)}

        # Step 2: Sort the second list based on the index in the first list
        sorted_genes = sorted(genes, key=lambda x: index_map.get(x, float('inf')))

        sorted_genes = [tf] + sorted_genes   

        # Create a PdfPages object
        # Create a PdfPages object
        if len(sorted_genes) <= 8:
            size = 6
        elif len(sorted_genes) > 8 and len(sorted_genes) <= 10:
            size = 8
        else:
            size = 10

        fig, axs = plt.subplots(1, 2, figsize=(10, size))
        sc.set_figure_params(scanpy=True, fontsize=20)
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


# # 6ho_S2

# In[31]:


#Read the stats dataframe
s2_6ho_stats = pd.read_csv(f'{stats_dir}6Ho/S2_6Ho.csv')
s2_6ho_stats   


# In[32]:


curr_enrichment_df = pd.read_csv(f'{code_dir}/Enrichment_tables/Stem_Cell_6ho_S2.tsv',sep='\t',index_col=0)

#for sort in sorts:
    
pdf_filename = os.path.join(save_dir, f"Stem_Cell_reg_6ho_S2_dp_sort_stats_v4.pdf")
with PdfPages(pdf_filename) as pdf:
    curr_metadata = metadata[metadata['cluster1'] == 'Mammary epithelial cells-Basal']
    curr_sub_info = [clust.split('l_')[1] for clust in curr_metadata['sub_info']]
    curr_metadata['curr_sub_info'] = curr_sub_info
    adata_scaled_basal.obs = curr_metadata

    # Filter for 5ht and 6ho
    adata_5ht_basal = adata_scaled_basal[adata_scaled_basal.obs['orig.ident'] == '5Ht']
    adata_6ho_basal = adata_scaled_basal[adata_scaled_basal.obs['orig.ident'] == '6Ho']


    adata_s1 = adata_6ho_basal[adata_6ho_basal.obs['curr_sub_info']=='S1']
    adata_s2 = adata_6ho_basal[adata_6ho_basal.obs['curr_sub_info']=='S2']
    adata_s3 = adata_6ho_basal[adata_6ho_basal.obs['curr_sub_info']=='S3']
    adata_s4 = adata_6ho_basal[adata_6ho_basal.obs['curr_sub_info']=='S4']
    adata_u1 = adata_6ho_basal[adata_6ho_basal.obs['curr_sub_info']=='U1_ko']
    adatas = [adata_s1,adata_s2,adata_s3,adata_s4,adata_u1]

    for i in range(len(curr_enrichment_df)):
        regulon = curr_enrichment_df['Transcription Factor'][i]
        tf = regulon.split('(')[0]
        genes = curr_enrichment_df['DETG'][i].split(', ')
        genes = [gene for gene in genes if gene != '']
        # Ensure the directory exists
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)


            # Sort by order from stats file

        # Extract the specific cell (which is a string representation of a list)
        target_genes_string = s2_6ho_stats[s2_6ho_stats['Significant_Regulons'] == regulon]['Significant Target Genes (5Ht_6Ho)'].iloc[0]

        # Convert the string to a list using ast.literal_eval
        target_genes_list = ast.literal_eval(target_genes_string)

        # Step 1: Create a mapping from the first list to its index position
        index_map = {value: index for index, value in enumerate(target_genes_list)}

        # Step 2: Sort the second list based on the index in the first list
        sorted_genes = sorted(genes, key=lambda x: index_map.get(x, float('inf')))

        sorted_genes = [tf] + sorted_genes   

        # Create a PdfPages object
        if len(sorted_genes) <= 8:
            size = 6
        elif len(sorted_genes) > 8 and len(sorted_genes) <= 10:
            size = 8
        else:
            size = 10

        # Create a PdfPages object
        fig, axs = plt.subplots(1, 2, figsize=(10, size))
        sc.set_figure_params(scanpy=True, fontsize=20)
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

# In[5]:

# single out 'Hif1a', 'Stat3' 
# #Always keep the genes to be single out at the top in the desired order
s2_stem_cell_diff = ['Hif1a', 'Stat3', 'Fgfr1', 'Rps7', 'Msi2', 'Mef2c', 'Nrg1', 'Sema6d', 'Sox9', 'Epcam', 
                     'Kitl', 'Zfp36l2', 'Anxa6',  'Zfp36', 'Nfe2l2', 'Lrp6', 'Rbbp4', 'Mtch2', 'Rbbp7']


# single out 'Stat3'
s2_stem_cell_pop = [
    'Stat3','Smarca2', 'Fgfr1', 'Arid4b', 'Jmjd1c', 'Sox9', 'Elavl1',
    'Zfp36l2', 'Crebbp', 'Ss18', 'Ddx6', 'Smc1a', 'Rbbp4', 'Rbbp7'
]


s2_stem_cell_prol = [
    'Fgfr1', 'Zfp36l1', 'Etv6', 'Hmgb2', 'Sox9', 'Epcam', 'Kitl',
    'Pbx1', 'Lrp6', 'Cxcl1'
]

u1_stem_cell_diff = [
    'Ptprc', 'Zeb2', 'Fn1', 'Nrp2', 'Ap2a2', 'Jarid2', 'Nrp1', 'Chd2',
    'Cited2', 'Coro1c', 'Zfp36', 'Sema4a', 'Cdc42', 'Mapk14', 'Rbpj',
    'Ncoa3', 'Vsir', 'Mef2c'
]

all_s2_stem_cell = [s2_stem_cell_diff,s2_stem_cell_pop,s2_stem_cell_prol]

single_out_num = 2

# # U1_wt

# In[7]:


sorts = ['Fraction','Mean']

for sort in sorts:
    pdf_filename = os.path.join(save_dir, f"Stem_Cell_DEG_U1_wt_dp_sort_{sort}_v4.pdf")
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


        #genes = u1_stem_cell_diff
        # Ensure the directory exists
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        if sort == 'Fraction':

            fraction_expressed = (adata_u1[:, u1_stem_cell_diff].X > 0).mean(axis=0)

            # Convert to a pandas Series for sorting
            fraction_expressed_series = pd.Series(fraction_expressed, index=u1_stem_cell_diff)

            # Step 2: Sort genes by fraction of cells
            sorted_genes = fraction_expressed_series.sort_values(ascending=False).index.tolist()

        else:
            gene_exp = []
            for g in u1_stem_cell_diff:
                mean_exp_gene = []
                for data in adatas:
                    mean_exp_gene.append(data[:, [g]].to_df().mean(axis=0)[0])
                mean_exp_gene = (mean_exp_gene - np.min(mean_exp_gene)) / (np.max(mean_exp_gene) - np.min(mean_exp_gene))
                gene_exp.append(mean_exp_gene[-1])#the -1 index represents cluster u1_wt

            # Step 1: Pair the elements together
            paired_lists = list(zip(u1_stem_cell_diff, gene_exp))

            # Step 2: Sort based on the second list (list2), in descending order
            sorted_pairs = sorted(paired_lists, key=lambda x: x[1], reverse=True)

            # Step 3: Unzip the sorted pairs back into two lists
            sorted_list1, sorted_list2 = zip(*sorted_pairs)

            # Convert them back to lists (since zip returns tuples)
            sorted_genes = list(sorted_list1)


 
        # Create a PdfPages object
        if len(sorted_genes) <= 8:
            size = 6
        elif len(sorted_genes) > 8 and len(sorted_genes) <= 10:
            size = 8
        else:
            size = 10

        # Create a PdfPages object
        fig, axs = plt.subplots(1, 2, figsize=(10, size))
        sc.set_figure_params(scanpy=True, fontsize=20)
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


# # S2

# In[ ]:


sorts = ['Fraction','Mean']
types = ['stem_cell_diff','stem_cell_maint','stem_cell_prolif']
for t,s2_DEGs in enumerate(all_s2_stem_cell):
    for sort in sorts:
        pdf_filename = os.path.join(save_dir, f"Stem_Cell_DEG_S2_dp_sort_{types[t]}_{sort}_v4.pdf")
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

            if s2_DEGs == s2_stem_cell_diff:#Single out Hif1 and Stat3 in stem cell diff pathway DEGs;sort the rest
                single_out_num = 2 #Change as per requirement. 
            elif s2_DEGs == s2_stem_cell_pop:#Single out Stat3 in stem cell pop maintainence pathway DEGs;sort the rest
                single_out_num = 1 #Change as per requirement.
            else:
                single_out_num = 0

            if single_out_num > 0:
                singled_out_genes = s2_DEGs[:single_out_num]
                remaining_genes = s2_DEGs[single_out_num:]
            else:
                remaining_genes = s2_DEGs

            #genes = u1_stem_cell_diff
            # Ensure the directory exists
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)

            if sort == 'Fraction':

                fraction_expressed = (adata_s2[:, remaining_genes].X > 0).mean(axis=0)

                # Convert to a pandas Series for sorting
                fraction_expressed_series = pd.Series(fraction_expressed, index=remaining_genes)

                # Step 2: Sort genes by fraction of cells
                sorted_genes = fraction_expressed_series.sort_values(ascending=False).index.tolist()

            else:
                gene_exp = []
                for g in remaining_genes:
                    mean_exp_gene = []
                    for data in adatas:
                        mean_exp_gene.append(data[:, [g]].to_df().mean(axis=0)[0])
                    mean_exp_gene = (mean_exp_gene - np.min(mean_exp_gene)) / (np.max(mean_exp_gene) - np.min(mean_exp_gene))
                    gene_exp.append(mean_exp_gene[1])#the 1 index represents cluster S2

                # Step 1: Pair the elements together
                paired_lists = list(zip(remaining_genes, gene_exp))

                # Step 2: Sort based on the second list (list2), in descending order
                sorted_pairs = sorted(paired_lists, key=lambda x: x[1], reverse=True)

                # Step 3: Unzip the sorted pairs back into two lists
                sorted_list1, sorted_list2 = zip(*sorted_pairs)

                # Convert them back to lists (since zip returns tuples)
                sorted_genes = list(sorted_list1)

            if single_out_num > 0:
                sorted_genes = singled_out_genes + sorted_genes

            # Create a PdfPages object
            if len(sorted_genes) <= 8:
                size = 6
            elif len(sorted_genes) > 8 and len(sorted_genes) <= 10:
                size = 8
            else:
                size = 10
            # Create a PdfPages object
            fig, axs = plt.subplots(1, 2, figsize=(10, size))
            
            sc.set_figure_params(scanpy=True, fontsize=20)
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


