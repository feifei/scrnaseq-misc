import requests
import pandas as pd
from io import StringIO


def request_post(url, identifiers, species=10090, caller_identity='feifei_uu'):
    """ post requests to API """
    params = {'identifiers': identifiers,
              'species': species,
              'caller_identity': caller_identity}
    r = requests.post(url, data=params)
    return r
    
    
def stringdb(genes_list, db_version='11-0b', species=10090, caller_identity='feifei_uu'):
    ''' Run stringdb and return enrichment analysis table, 
        Newest version 11-5, we were using 11-0b,
        Default species is mouse
    '''
     
    url = 'https://version-%s.string-db.org' %db_version #'https://string-db.org' 
    get_string_ids_url = url + '/api/tsv/get_string_ids?'
    enrichment_url = url + '/api/tsv/enrichment?'
    
    gene_identifiers = '%0d'.join(genes_list)
    r = request_post(get_string_ids_url, gene_identifiers, species=species, caller_identity=caller_identity)
    df_stingid = pd.read_csv(StringIO(r.text), sep='\t')
    stringid_identifiers = '%0d'.join(list(df_stingid['stringId']))
    r = request_post(enrichment_url, stringid_identifiers, species=species, caller_identity=caller_identity)
    df = pd.read_csv(StringIO(r.text), sep='\t')
    return df
    

def adjust_columns(df):
    """ Adjust columns so that it's inline with string-db output """
    new_df = pd.DataFrame()
    new_df['#term ID'] = df['term']
    new_df['term description'] = df['description']
    new_df['observed gene count'] = df['number_of_genes']
    new_df['background gene count'] = df['number_of_genes_in_background']
    new_df['p_value'] = df['p_value']
    new_df['false discovery rate'] = df['fdr']
    new_df['matching proteins in your network (IDs)'] = df['inputGenes']
    new_df['matching proteins in your network (labels)'] = df['preferredNames']
    return new_df


def save_df(df, dir, prefix, tag = ''):
    """ save enrichment output in desired format, and sub files """
    
    outfile = dir + prefix + '.string%s.csv' %tag
    # all
    df.to_csv(outfile, index=False, float_format='%.2e')
    # bp
    sub_df = df[df['category'] == 'Process']
    sub_df = adjust_columns(sub_df)
    outfile = dir + prefix + '.string%s.bp.csv' %tag
    sub_df.to_csv(outfile, index=False, float_format='%.2e')
    # cc
    sub_df = df[df['category'] == 'Component']
    sub_df = adjust_columns(sub_df)
    outfile = dir + prefix + '.string%s.cc.csv' %tag
    sub_df.to_csv(outfile, index=False, float_format='%.2e')
    # bp
    sub_df = df[df['category'] == 'Function']
    sub_df = adjust_columns(sub_df)
    outfile = dir + prefix + '.string%s.mf.csv' %tag
    sub_df.to_csv(outfile, index=False, float_format='%.2e')
    # bp
    sub_df = df[df['category'] == 'KEGG']
    sub_df = adjust_columns(sub_df)
    outfile = dir + prefix + '.string%s.kegg.csv' %tag
    sub_df.to_csv(outfile, index=False, float_format='%.2e')