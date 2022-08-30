import scanpy as sc
import gseapy as gp
import matplotlib.pyplot as plt
from plot_func import *
from run_string import *

def rank_genes(adata, qry, ref, target_group, use_raw=True, method='wilcoxon', n_genes=50):
    ''' Rank genes based on qry, ref, target_group '''
    # DE should be calculated on all genes instead of HVG
    # Also on normalized and logrithimized values, instead of scaled values
    sc.tl.rank_genes_groups(adata, target_group, groups=[qry], 
                            reference=ref, method=method, use_raw=use_raw) 
    with plt.rc_context({"figure.figsize": (10, 5)}):
        sc.pl.rank_genes_groups(adata, groups=[qry], n_genes=n_genes)
    #sc.pl.rank_genes_groups_violin(adata, groups=qry, n_genes=n_genes) #use_raw=True by default
    # Plot the same genes as violins across all datasets.



def enr_analysis(glist, gene_sets, description, outfile, cutoff=0.05):
    ''' Enrichment analysis'''
    enr = gp.enrichr(gene_list=glist, organism='Mouse',
                     gene_sets=gene_sets,
                     description=description, cutoff=cutoff)
    if gene_sets == 'GO_Biological_Process_2018':
        with plt.rc_context({"figure.figsize": (4, 4)}):
            gp.plot.barplot(enr.res2d, title=description, cutoff=cutoff)
    results = enr.results[enr.results['P-value'] <= cutoff]
    if results.shape[0] >0:
        results.to_csv(outfile, index=False)


def multi_enr_analysis(glist, work_dir, qry, ref, extra_tag, cutoff, up_down=''):
    ''' Call enr_analysis, write output '''
    desc = '%s_vs_%s.%skegg%s' %(qry, ref, extra_tag, up_down)
    outfile = '%s%s.csv' %(work_dir, desc)
    enr_analysis(glist, 'KEGG_2019_Mouse', desc, outfile, cutoff)
    desc = '%s_vs_%s.%sgo_mf%s' %(qry, ref, extra_tag, up_down)
    outfile = '%s%s.csv' %(work_dir, desc)
    enr_analysis(glist, 'GO_Molecular_Function_2018', desc, outfile, cutoff)
    desc = '%s_vs_%s.%sgo_cc%s' %(qry, ref, extra_tag, up_down)
    outfile = '%s%s.csv' %(work_dir, desc)
    enr_analysis(glist, 'GO_Cellular_Component_2018', desc, outfile, cutoff)
    desc = '%s_vs_%s.%sgo_bp%s' %(qry, ref, extra_tag, up_down)
    outfile = '%s%s.csv' %(work_dir, desc)
    enr_analysis(glist, 'GO_Biological_Process_2018', desc, outfile, cutoff)


        
def get_glist(sig_genes):
    ''' Get gene list from sig_genes '''
    if sig_genes.shape[0] >1:
        glist = sig_genes['names'].squeeze().str.strip().tolist()
    elif sig_genes.shape[0] == 1:
        glist = [sig_genes['names'].squeeze().strip()]
    else:
        glist = []
    
    return glist



def rank_and_enr(adata, qry, ref, target_group, extra_tag='', save_volcano=False, use_raw=True, work_dir='./data/', method='wilcoxon', n_genes=50, cutoff=0.05):
    ''' Wrapper to rank genes and do enrichment analysis 
        extra_tag to specify cell_type. and other information about the comparison
    '''
#    if use_raw:
#        # Default behavior, do not want extra tag
#        extra_tag = extra_tag + ''
#    else:
#        # this will use the hvg scaled data
#        extra_tag = extra_tag + 'hvg.'
    extra_tag=''
    work_dir_csv = work_dir + 'results/csv/'
    rank_genes(adata, qry, ref, target_group, use_raw, method, n_genes)
    all_genes = sc.get.rank_genes_groups_df(adata, group=qry)
    #all_genes = all_genes[all_genes['names'].isin(proteins['gene_name'].tolist())]
    all_genes = all_genes[~all_genes['names'].str.startswith('Gm')]
    sig_genes = all_genes[all_genes['pvals_adj'] <= cutoff]
    sig_genes.sort_values(by=['logfoldchanges'], ascending=False, inplace=True)
    sig_genes.to_csv('%s%s_vs_%s.%ssig_genes.csv' %(work_dir_csv, qry, ref, extra_tag), index=False)
    
    
    # dotplots and stackplot of top 50 and bottom 50
    up50 = sig_genes.head(n_genes)['names'].tolist()
    down50 = sig_genes.tail(50)['names'].tolist()
    down50.reverse()
    sc.pl.dotplot(adata, up50, groupby=target_group, title='%s vs. %s up' %(qry, ref))
    #sc.pl.tracksplot(adata, up50, groupby=target_group)
    sc.pl.dotplot(adata, down50, groupby=target_group, title='%s vs. %s down' %(qry, ref))
    #sc.pl.tracksplot(adata, down50, groupby=target_group)



    glist = get_glist(sig_genes)
    # Separate up and downs
    sig_genes_up = sig_genes[sig_genes['logfoldchanges'] > 0]
    sig_genes_down = sig_genes[sig_genes['logfoldchanges'] < 0]
    glist_up = get_glist(sig_genes_up)
    glist_down = get_glist(sig_genes_down)

    print('%d sig_genes, with %d up and %d down' %(len(glist), len(glist_up), len(glist_down)))
    
    if len(glist) == 0:
        return None
        #multi_enr_analysis(glist, work_dir, qry, ref, extra_tag, desc)
    
    # volcano plot
    # plt.rcParams["figure.figsize"] = [10,5]
    # sig_genes.plot.scatter(x='logfoldchanges', y='pvals_adj', logy=True)
    # plt.gca().invert_yaxis()
    top_up_genes = sig_genes.head(20)['names'].tolist()
    top_down_genes = sig_genes.tail(20)['names'].tolist()
    plot_volcano(all_genes, top_up_genes+top_down_genes, adjust=True, figsize=(5, 8))
    
    if save_volcano:
        work_dir_plot = work_dir + 'results/plot/'
        volcano_file = '%s/%s_vs_%s.volcano.pdf' %(work_dir_plot, qry, ref)
        plt.savefig(volcano_file, dpi=300)
    
    #plt.rcParams["figure.figsize"] = [4,4] 
    
    if len(glist_up) > 0:
        multi_enr_analysis(glist_up, work_dir_csv, qry, ref, extra_tag, cutoff, up_down='-up')
    if len(glist_down) > 0:
        multi_enr_analysis(glist_down, work_dir_csv, qry, ref, extra_tag, cutoff, up_down='-down')

        
    return sig_genes
       



def rank_and_string(adata, qry, ref, target_group, proteins, extra_tag='', save_volcano=False, use_raw=True, 
                    work_dir='./data/', method='wilcoxon', n_genes=50, cutoff=0.05):
    ''' Wrapper to rank genes and do enrichment analysis with string
        extra_tag to specify cell_type. and other information about the comparison
    '''
    if use_raw:
        # Default behavior, use raw data
        extra_tag = extra_tag + ''
    else:
        # this will use the hvg scaled data
        extra_tag = extra_tag + 'hvg.'
    #work_dir_csv = work_dir + 'results/csv_string/'
    prefix = '%s_vs_%s%s' %(qry, ref, extra_tag)
    rank_genes(adata, qry, ref, target_group, use_raw, method, n_genes)
    all_genes = sc.get.rank_genes_groups_df(adata, group=qry)
    # keep only protein coding genes, implemented for large intestine project due to rRNAs
    all_genes = all_genes[all_genes['names'].isin(proteins['gene_name'].tolist())]
    all_genes = all_genes[~all_genes['names'].str.startswith('Gm')]
    sig_genes = all_genes[all_genes['pvals_adj'] <= cutoff]
    #sig_genes.sort_values(by=['logfoldchanges'], ascending=False, inplace=True)
    sig_genes.to_csv('%s/%s.sig_genes.csv' %(work_dir, prefix), index=False)
    
    glist = get_glist(sig_genes)
    # Separate up and downs
    sig_genes_up = sig_genes[sig_genes['logfoldchanges'] > 0]
    sig_genes_down = sig_genes[sig_genes['logfoldchanges'] < 0]
    glist_up = get_glist(sig_genes_up)
    glist_down = get_glist(sig_genes_down)

    print('%d sig_genes, with %d up and %d down' %(len(glist), len(glist_up), len(glist_down)))
    
    if len(glist) == 0:
        return None
 
     
    # dotplots and stackplot of top 50 and bottom 50
    up50 = sig_genes.head(n_genes)['names'].tolist()
    down50 = sig_genes.tail(n_genes)['names'].tolist()
    down50.reverse()
    with plt.rc_context({"figure.figsize": (10, 5)}):
        sc.pl.dotplot(adata, up50, groupby=target_group, title='%s vs. %s up' %(qry, ref))
        #sc.pl.tracksplot(adata, up50, groupby=target_group)
        sc.pl.dotplot(adata, down50, groupby=target_group, title='%s vs. %s down' %(qry, ref))
        #sc.pl.tracksplot(adata, down50, groupby=target_group)



   
    # volcano plot, name the top 20 up and down genes
    top_up_genes = sig_genes.head(20)['names'].tolist()
    top_down_genes = sig_genes.tail(20)['names'].tolist()
    plot_volcano(all_genes, top_up_genes+top_down_genes, adjust=True, figsize=(5, 8))
    
    if save_volcano:
        work_dir_plot = './figures/'
        volcano_file = '%s/%s.volcano.pdf' %(work_dir_plot, prefix)
        plt.savefig(volcano_file, dpi=300)
    

    # Go enrichment analysis with string
    if len(glist_up) > 0:
        df_up = stringdb(glist_up)
        save_df(df_up, work_dir, prefix, tag = '-up')
    if len(glist_down) > 0:
        df_down = stringdb(glist_down)
        save_df(df_down, work_dir, prefix, tag = '-down')
        
    return sig_genes

