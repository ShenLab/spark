'''
16/08/2019 @queenjobo @ksamocha

Script to run enrichment test part of the DeNovoWEST testing framework

For usage see 'python DNE_test.py -h'

make a change
'''

__version__ = 1.2
__date__ = '2019-01-09'


#IMPORTS--------------

import logging
import random
import argparse
import pandas as pd
import numpy as np
import itertools
import datetime
from scipy import stats
from probabilities import *


#set seed
#random.seed(31)
#np.random.seed(31)
random.seed(2014)
np.random.seed(2014)

#get todays date
now = datetime.datetime.now()
today = str(now).split(' ')[0]

#default paths
wdicpath = "input/loess_weight_dic_v4.tab"
denovopath = "input/de_novos.txt.gz"
outpath = "DeNovoWEST_enrichment_results_" + today + ".tab"


#FUNCTIONS---------------

def get_options():
    '''
    parse options
    '''
    parser = argparse.ArgumentParser()
    #parser.add_argument('--weightdic', default=wdicpath,
    #    help='path to dictionary of weights to assign variants')
    parser.add_argument('--rates', help='path to all possible exome SNVs with specific rates and weights')
    parser.add_argument('--denovos', default=denovopath,  help='path to annotated denovo variants')
    parser.add_argument('--output', default= outpath, help = "output file to save to")
    parser.add_argument('--nmales',type = int, help = 'number of males in study')
    parser.add_argument('--nfemales', type = int, help = 'number of females in study')
    parser.add_argument('--hgbuild', default = "hg19", help = 'genome assembly build')
    parser.add_argument('--nsim',default = 1000000000, type = int,help = 'minimum number of simulations for each gene (default 1 billion)')
    parser.add_argument('--pvalcap',type = float, default = 1.0, help = "stop simulations if cumulative p-value > pvalcap (default 1)") 
    return parser.parse_args()

    
def load_denovos(path):
    ''' 
    get a DataFrame of DDD de novos
    '''
    logging.info("Loading de novos")
    denovos = pd.read_table(path, sep='\t', na_values='.', keep_default_na=True) # added sep
    denovos['Chrom'] = denovos['Chrom'].astype(str)
    denovos['Position'] = denovos['Position'].astype(int)
    # denovos['score'] = denovos['score'].astype(float)
    # denovos = fix_consequences(denovos)
    return(denovos)


def load_rates(ratespath,genes):
    '''
    load rates file
    '''
    #rates = pd.read_table(ratespath,names = ['symbol','chrom','pos','ref','alt','cq','prob','raw','score','maf','constrained'])
    logging.info("Loading rates")
    rates = pd.read_table(ratespath, sep='\t', na_values='.', keep_default_na=True)
    #subset to only genes observed in de novo file
    rates = rates[rates['GeneID'].isin(genes)] 
    rates['Chrom'] = rates['Chrom'].astype(str)
    #rates['Position'] = rates['Position'].astype(int)
    return(rates)

    
def correct_for_x_chrom(rates, male_n, female_n, hgbuild):
    ''' 
    correct for X chromosome
    '''
    autosomal = 2 * (male_n + female_n)
    female_transmit = male_n + female_n
    male_transmit = female_n
    
    # get scaling factors using the alpha from the most recent SFHS (Scottish
    # Family Health Study) phased de novo data.
    alpha = 3.4
    male_k = 2 / (1 + (1 / alpha))
    female_k = 2 / (1 + alpha)
    
    # correct the non-PAR chrX genes for fewer transmissions and lower rate
    # (dependent on alpha)
    # need to account for PAR part of chrX! <-- TBD
    if hgbuild=="hg19" or hgbuild=="b37":
        hapX = rates['Chrom'].isin(['X', 'chrX']) & rates['Position'].between(2699520,154931044)
    elif hgbuild=="hg38" or hgbuild=="b38":
        hapX = rates['Chrom'].isin(['X', 'chrX']) & rates['Position'].between(2781479,155701382)
    else:
        hapX = rates['Chrom'].isin(['X', 'chrX'])


    x_factor = ((male_transmit * male_k) + (female_transmit * female_k)) / autosomal
    x_factor = pd.Series([x_factor] * len(hapX), index=rates.index)
    x_factor[~hapX] = 1
    
    rates['prob'] *= x_factor
    
    return rates

    
def get_expected_rates(rates, male, female, hgbuild):
    '''
    get expected mutation rates 
    '''
    logging.info("get expected")
    autosomal = 2 * (male + female)
    rates['prob'] = rates['prob'] * autosomal
    
    return correct_for_x_chrom(rates, male, female, hgbuild)

    
#def test_gene(rates,denovos,gene,nsim,indel_weights,pvalcap):
def test_gene(rates,denovos,gene,nsim,pvalcap):
    '''
    test single gene given rates object and de novos
    '''
    generates = rates[rates.GeneID == gene]
    # Indel rates has been added to the rates df
    # generates = get_indel_rates(generates,indel_weights)
    # here add inframe/frameshift rates
    s_obs = denovos[denovos.GeneID == gene].weight.sum()
    pval,info,s_exp = get_pvalue(generates,s_obs,nsim,pvalcap)
    return(s_exp,s_obs,pval,info)


    
#def test_all_genes(rates,denovos,nsim,indel_weights,pvalcap):
def test_all_genes(rates,denovos,nsim,pvalcap):
    '''
    test all genes for enrichment
    '''
    logging.info("Starting tests: ")
    #genes = rates.symbol.unique()
    genes = denovos.GeneID.unique()
    results = []
    for gene in genes:
        #skip gene if no observed de novos in our dataset
        #if denovos[denovos.symbol == gene].weight.sum() == 0:
        #    continue
        #skip gene if no rates in our dataset
        if gene not in rates.GeneID.unique():
            logging.info("could not find " + str(gene))
            continue
        logging.info("testing " + str(gene))
        #s_exp,s_obs,pval,info = test_gene(rates,denovos,gene,nsim,indel_weights,pvalcap)
        s_exp,s_obs,pval,info = test_gene(rates,denovos,gene,nsim,pvalcap)
        # extract hgnc_id
        # hgnc = denovos.loc[denovos['GeneID'] == gene, 'hgnc_id'].iloc[0]
        # results.append((gene,hgnc,s_exp,s_obs,pval,info))
        results.append((gene,s_exp,s_obs,pval,info))
    return(results)

    
def main():
    args = get_options()

    #initialise log file
    logging.basicConfig(filename=args.output.split(".")[0]+".log",level=logging.DEBUG)

    #load denovos
    denovos = load_denovos(args.denovos)
    #get genes in de novo file so only run test on these
    genes = denovos.GeneID.unique()
    logging.info("genes:" + ",".join(genes))
    
    #use de novo mutation weights
    #denovos,_ = get_weights(denovos,args.weightdic,False)
    denovos = denovos[denovos.weight.apply(np.isreal)]
    denovos = denovos.dropna(subset=['weight'])
    
    #load rates and get expected
    rates = load_rates(args.rates,genes)
    rates = rates.dropna(subset=['weight','prob'])

    logging.info("rates loaded")
    if len(rates.index) == 0:
        logging.info("rates are empty")
        return

    if args.hgbuild != "hg19" and args.hgbuild != "b37" and args.hgbuild != "hg38" and args.hgbuild != "b38":
        logging.info("cannot recognize genome assembly build : " + args.hgbuild)
        return

    rates = get_expected_rates(rates, args.nmales, args.nfemales, args.hgbuild)
    #get weights
    #logging.info("get rates weights")
    #rates,indel_weights = get_weights(rates,args.weightdic,False)
    rates = rates[rates.weight.apply(np.isreal)]
    
    # run gene tests
    #results = test_all_genes(rates,denovos,args.nsim,indel_weights,args.pvalcap)
    results = test_all_genes(rates,denovos,args.nsim,args.pvalcap)

    #save results to file
    # df = pd.DataFrame.from_records(results,columns = ["symbol","hgnc_id","expected","observed","p-value","info"])
    df = pd.DataFrame.from_records(results,columns = ["GeneID","Expected","Observed","P-value","Info"])
    df.to_csv(args.output, sep = "\t", index = False)

#------------------SCRIPT--------------------------------------

if __name__=='__main__':
    main()
    
    
    
