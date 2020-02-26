import numpy as np
from scipy.stats import chi2,norm
import pandas as pd
import re
import argparse
import sys
from cleangwas import read_header,parse_header,default_cnames, timeit,write_file,get_converter
from decimal import Decimal,getcontext
import logging
import os
import mpmath
from platform import system
from subprocess import getstatusoutput
from biomart import BiomartServer


'''
1. p值 z值互相转换
   标记Z值方向
   beta se 转 P
   beta se 转 Z
   分列
2. OR值 BETA值相互转换。对应的SE也会转换。TODO
3. 更新chromesome position 
4. 根据CHR POS 更新RS
5. 添加frequency
'''

mpmath.mp.prec = 3000
getcontext().prec = 20
class options():
    def __init__(self):
        self.infp = 'test.gwas'
        self.outfp = 'test'
        self.sep = '\s+'
        self.beta_to_p = True
        self.split=None
        self.beta_to_z=False
        self.p_to_z = 'sign'
        self.z_to_p = False
        self.convert_coordinate = '37to38'
        self.rs_name = None
        self.chr_name = None
        self.pos_name = None
        self.a1_name = None
        self.a2_name = 'Aa2'
        self.freq_name = 'FREQ.A1.1000G.EUR'
        self.beta_name = None
        self.or_name = None
        self.se_name = None
        self.info_name = None
        self.P_name = None
        self.N_name = None
class Logger():
    def __init__(self, name=None):
        self.logger = logging.getLogger(name) if name else logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        self.formatter = logging.Formatter('%(message)s')
        self.f_handler = logging.FileHandler('{0}.log'.format(__name__))
        self.f_handler.setFormatter(self.formatter)
        self.f_handler.setLevel(logging.INFO)
        self.logger.addHandler(self.f_handler)
        self.console_handler = logging.StreamHandler()
        self.console_handler.setFormatter(self.formatter)
        self.console_handler.setLevel(logging.INFO)
        self.logger.addHandler(self.console_handler)
        self.logger.setLevel(logging.INFO)

    def rename(self, outprefix):
        self.f_handler.flush()
        self.f_handler.close()
        self.logger.removeHandler(self.f_handler)
        if isinstance(outprefix, str):
            try:
                os.remove('{0}.log'.format(outprefix))
            except:
                pass
            try:
                os.rename('{0}.log'.format(__name__), '{0}.log'.format(outprefix))
            except:
                print('Can\'t rename log file: {0}'.format('{0}.log'.format(__name__)))
logger = Logger()


def parse_arg():
    parser = argparse.ArgumentParser(description='Transform GWAS data.')
    parser.add_argument('-i', '--input', dest='infp', required=True, default=None, metavar='FILE', help='Input file')
    parser.add_argument('-o', '--output', dest='outfp', required=True, default=None, metavar='FILE', help='Output file')
    parser.add_argument('--sep', default='\s+', metavar='separator', help='delimiter of input file')
    parser.add_argument('--split', default=None, nargs=2, action='append',
                        help='Split columns. e.g., SNPID "[chr@CHR]:@[POS]"')
    parser.add_argument('--beta-to-p', default=None, choices=['oneside','twoside'], help='Calculate P value from beta and se.')
    parser.add_argument('--beta-to-z', default=False, action='store_true', help='Calculate Z value from beta and se.')
    parser.add_argument('--p-to-z', default=None, choices=['sign', 'unsign'], help='Calculate Z value from P value. Z will be signed if --p-to-z was specified with sign')
    parser.add_argument('--z-to-p', default=False, action='store_true',  help='Calculate P value from Z value')
    parser.add_argument('--convert-coordinate', default=None, choices=['36to37', '36to38', '37to36', '37to38', '38to37'],
                        help='Convert coordinates between different assemblies.')

    parser.add_argument('--rs-name', help='Specify rs column name')
    parser.add_argument('--chr-name', help='Specify chromosome column name')
    parser.add_argument('--pos-name', help='Specify position column name')
    parser.add_argument('--a1-name', help='Specify effect allele column name')
    parser.add_argument('--a2-name', help='Specify non effect allele column name')
    parser.add_argument('--freq-name', help='Specify frequency column name')
    parser.add_argument('--beta-name', help='Specify beta coefficient column name')
    parser.add_argument('--or-name', help='Specify OR coefficient column name')
    parser.add_argument('--se-name', help='Specify SD column name')
    parser.add_argument('--info-name', help='Specify imputation score column name')
    parser.add_argument('--N-name', help='Specify sample size column name')
    parser.add_argument('--P-name', help='Specify P value column name')
    return parser.parse_args()


def p_to_z(p):
    return np.sqrt(chi2.isf(p, 1))

def signz(z, sign, value=0):
    df = pd.DataFrame(zip(abs(z), sign), columns=['Z', 'SIGN'])
    return df.apply(lambda x: x['Z'] if x['SIGN'] > value else -1 * x['Z'], axis=1)

def z_to_p(z, oneside=False):
    tail = 1 if oneside else 2
    return ((-abs(z)).apply(mpmath.ncdf) * tail).apply(mpmath.nstr, n=6)

def beta_to_p(beta, se, oneside=False):
    return z_to_p(beta/se, oneside)

def beta_to_z(beta, se):
    return beta/se

def split_col(df, name, ex, keep=False):
    '''
    ex expample:  [@RS]:[chr@CHR]:[@POS]
    '''
    fields = re.findall('@\w+]', ex)
    fieldnames = [x.strip('@').strip(']') for x in fields]
    newex = ex.replace('*', '\*')
    for x in fields: newex = newex.replace(x, '(.*)').replace('[', '')
    expr = re.compile(newex, re.I)
    newseries = df[name].map(lambda x: '*'.join(re.match(expr, x).groups()))
    expanddf = newseries.str.split(pat='*', expand=True)
    expanddf.columns = fieldnames
    if not keep:
        df.drop(columns=name, inplace=True, axis=1)
    df = pd.concat([df, expanddf], axis=1,  sort=False)
    return df

def liftover(df, conv_flag, chr, pos):
    pf = system()
    if not pf == 'Linux':
        print('This function can not be executed in {0}. '
              'Please use it in Linux platform\nExit...'.format(pf))
        sys.exit(1)
    if set([chr, pos]) - set(df.columns):
        logger.logger.warning('\033[5;31mWarning: genomic cooridnates can not be updated,'
                              ' need {0} columns. data will be not changed.\033[0m'.format(set([chr, pos]) - set(df.columns)))
        return df
    newdf = df.astype(str)
    chain_dict ={
        '36to37': 'hg18ToHg19.over.chain.gz',
        '36to38': 'hg18ToHg38.over.chain.gz',
        '37to36': 'hg19ToHg18.over.chain.gz',
        '37to38': 'hg19ToHg38.over.chain.gz',
        '38to37': 'hg38ToHg19.over.chain.gz',
    }
    newdf[chr] = newdf[chr].apply(lambda x : re.sub(r'^chr', 'chr', x, flags=re.I))
    cols = [chr, pos]
    remaincols = [x for x in list(newdf.columns) if not x in cols]
    newdf[','.join(remaincols)] = newdf.apply(lambda x :','.join(x[remaincols]), axis=1)
    newdf = pd.concat([df[[chr,pos]], newdf[','.join(remaincols)]], axis=1)
    newdf.insert(2, '__POS2__', df[pos]+1)
    newdf.to_csv('liftover.tmp', sep='\t', na_rep='NA',index=False, header=False)
    code,out = getstatusoutput('{0} liftover.tmp {1} .map .unmap'.format(
        os.path.join(os.getcwd(), 'refsource', 'liftover', 'liftOver'),
        os.path.join(os.getcwd(), 'refsource', 'liftover', chain_dict.get(conv_flag))))
    if code:
        logger.logger.warning('\033[5;31mWarning: Error in converting coordinates! coordinates of data will be not changed \033[0m')
        return df
    os.remove('liftover.tmp')
    usecols = list(newdf.columns)
    usecols.remove('__POS2__')
    newdf = pd.read_csv('.map', header=None, names=newdf.columns, dtype=str, usecols=usecols, sep='\t')
    split_name = newdf.columns[2].split(',')
    split_df = newdf.iloc[:, 2].str.split(',', expand=True)
    split_df.columns = split_name
    newdf = pd.concat([newdf.loc[:,[chr, pos]], split_df], axis=1)
    newdf[chr] = newdf[chr].apply(lambda x : re.sub(r'^chr', '', x, flags=re.I))
    '''统计失败SNP个数。
    with open('.unmap', 'w') as reader:
        lines = reader.readlines()
    delno = len({x for x in lines if not x.startswith('#')})
    '''
    os.remove('.map')
    os.remove('.unmap')
    return newdf

def calcu(df, func, overwrite:str, cnames:list, *use:tuple, **kargs):
    if overwrite in cnames:
        logger.logger.info('\033[5;31mWarning:\033[0m {0} column are already present in'
        'your file. new estimated {0} value will overwrite it.'.format(overwrite))
    if set(use) - set(cnames):
        logger.logger.warning('\033[5;31mWarning: {0} can not be estimated,'
                              ' need {1} columns. \033[0m'.format(overwrite, set(use)-set(cnames)))
        return df
    use_series = []
    for x in use: use_series.append(df[x])
    if len(kargs) > 0:
        df[overwrite] = func(*tuple(use_series), **kargs)
    else:
        df[overwrite] = func(*tuple(use_series))
    return df

def mainbody(opts, cnames, logger):
    total_df = pd.DataFrame()
    for chunk in pd.read_csv(opts.infp, sep=opts.sep, header=0, names=cnames, dtype=str, iterator=True,
                             chunksize=2000000):
        '''不在read_csv()使用converters参数是因为当文件最后一列是空值时，会报decimal.InvalidOperation: [<class 'decimal.ConversionSyntax'>]。
        另外同时用names和converters选项会报警告。'''
        for k, v in get_converter(cnames).items(): chunk[k] = chunk[k].apply(v)
        if opts.beta_to_p:
            chunk = calcu(chunk, beta_to_p, 'P', cnames, 'BETA', 'SE', oneside=opts.beta_to_p == 'oneside')
        if opts.beta_to_z:
            chunk = calcu(chunk, beta_to_z, 'Z', cnames, 'BETA', 'SE')
        if opts.p_to_z == 'sign':
            if 'BETA' in cnames:
                chunk = calcu(chunk, p_to_z, 'Z', cnames, 'P')
                chunk = calcu(chunk, signz, 'Z', cnames, 'Z', 'BETA', value=0)
            elif 'OR' in cnames:
                chunk = calcu(chunk, p_to_z, 'Z', cnames, 'P')
                chunk = calcu(chunk, signz, 'Z', cnames, 'Z', 'OR', value=1)
            else:
                logger.logger.warning('\033[5;31mWarning: Can not sign Z value in absence of beta or or statistic\033[0m \nExit...')
                sys.exit(1)
        elif opts.p_to_z == 'unsign':
            chunk = calcu(chunk, p_to_z, 'Z', cnames, 'P')
        if opts.z_to_p:
            chunk = calcu(chunk, z_to_p, 'P', cnames, 'Z')
        if opts.convert_coordinate:
            chunk = liftover(chunk, opts.convert_coordinate, 'CHR', 'POS')
        if opts.split:
            for x in opts.split:
                chunk = split_col(chunk, x[0], x[1], keep=True)
        total_df = total_df.append(chunk, sort=False)
    return total_df

@timeit(logger)
def init(opts, logger):
    header = read_header(opts.infp, opts.sep)
    cnames = parse_header(header, opts, logger, default_cnames)
    df = mainbody(opts, cnames,logger)
    write_file(df, opts, cnames)
    print(df)

if __name__ == '__main__':
    if len(sys.argv) == 1:
        opts = options()
    else:
        opts = parse_arg()
    init(opts, logger)
    logger.logger.info('End')
    logger.rename(opts.outfp)
    sys.exit(0)