#!/usr/bin/env python

#  This file is part of Tome
#
#  Tome is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Tome is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Tome.  If not, see <https://www.gnu.org/licenses/>


import sys
import os
import pandas as pd
from Bio import SeqIO
# 使用标准的 joblib 导入
import joblib
from collections import Counter
from multiprocessing import Pool, cpu_count
import numpy as np
import subprocess
import requests  # 添加 requests 用于跨平台下载
import shutil   # 添加 shutil 用于跨平台文件操作
import platform # 添加 platform 用于检测操作系统
import warnings

# 忽略一些兼容性警告
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)


# Esimtation of OGT for organism(s)
################################################################################
def print_out(line):
    sys.stdout.write(str(line)+'\n')

def parse_args():
    args = dict()
    for i in range(len(sys.argv)):
        item = sys.argv[i]
        if item.startswith('-'):
            item = item.replace('-','')
            try:args[item] = sys.argv[i+1]
            except:args[item] = ''

            if item == 'h': args['help'] = ''

    for i in range(len(sys.argv)):
        if 'tome' in sys.argv[i]:
            try:
                if sys.argv[i+1] in ['predOGT','getEnzymes']:
                    args['method'] = sys.argv[i+1]
            except: None
            break
    return args

def load_means_stds(predictor):
    means=dict()
    stds=dict()
    features=list()
    for line in open(predictor.replace('pkl','f'),'r'):
        if line.startswith('#'):continue
        cont=line.split()
        means[cont[0]]=float(cont[1])
        stds[cont[0]]=float(cont[2])
        features.append(cont[0])
    return means,stds,features

def train_model():
    from sklearn import svm
    from sklearn.metrics import r2_score
    from scipy.stats import spearmanr, pearsonr
    from sklearn.metrics import mean_squared_error
    from sklearn.preprocessing import StandardScaler

    path = os.path.dirname(os.path.realpath(__file__))
    predictor = os.path.join(path,'model/OGT_svr.pkl')
    
    def get_standardizer(X):
        # 使用 sklearn 的 StandardScaler 替代手动计算
        scaler = StandardScaler()
        scaler.fit(X)
        return scaler.mean_, scaler.scale_

    def standardize(X):
        # 使用 sklearn 的 StandardScaler 进行标准化
        scaler = StandardScaler()
        return scaler.fit_transform(X)

    # load training dataset
    trainfile = os.path.join(path,'data/train.csv')
    df = pd.read_csv(trainfile, index_col=0)
    X = df.values[:,:-1]
    Y = df.values[:,-1].ravel()
    features = df.columns[:-1]

    # 使用 StandardScaler 进行标准化
    scaler = StandardScaler()
    Xs = scaler.fit_transform(X)
    
    # 使用更新的 SVR 参数（scikit-learn 1.7.0 兼容）
    model = svm.SVR(kernel='rbf', C=64.0, epsilon=1.0, gamma='scale')
    model.fit(Xs, Y)

    # get model performance:
    p = model.predict(Xs)
    rmse = np.sqrt(mean_squared_error(Y, p))
    r2 = r2_score(Y, p)
    r_spearman = spearmanr(p, Y)
    r_pearson = pearsonr(p, Y)

    print_out('A new model has been successfully trained.')
    print_out('Model performance:')
    print_out('        RMSE: '+ str(rmse))
    print_out('          r2: ' + str(r2))
    print_out('  Pearson r: ' + str(r_pearson))
    print_out('  Spearman r: ' + str(r_spearman))
    print_out('')

    # save model and scaler
    print_out('Saving the new model to replace the original one...')
    joblib.dump(model, predictor)
    
    # 保存标准化参数
    fea = open(predictor.replace('pkl','f'),'w')
    means, stds = scaler.mean_, scaler.scale_
    fea.write('#Feature_name\tmean\tstd\n')
    for i in range(len(means)):
        fea.write('{0}\t{1}\t{2}\n'.format(features[i], means[i], stds[i]))
    fea.close()
    print_out('Done!')
    print_out('')

def load_model():
    path = os.path.dirname(os.path.realpath(__file__))
    predictor = os.path.join(path,'model/OGT_svr.pkl')
    try:
        # 使用 joblib 加载模型
        model = joblib.load(predictor)
        means, stds, features = load_means_stds(predictor)
    except Exception as e:
        print_out(f'Failed loading the model: {str(e)}. Trying to train the model...')
        train_model()
        model = joblib.load(predictor)
        means, stds, features = load_means_stds(predictor)

    return model, means, stds, features

def do_count(seq):
    dimers = Counter()
    for i in range(len(seq)-1): dimers[seq[i:i+2]] += 1.0
    return dimers


def count_dimer(fasta_file,p):
    seqs = [str(rec.seq).upper() for rec in SeqIO.parse(fasta_file,'fasta')]

    if p == 0:num_cpus = cpu_count()
    else: num_cpus = p
    results = Pool(num_cpus).map(do_count, seqs)
    dimers = sum(results, Counter())
    return dict(dimers)

def get_dimer_frequency(fasta_file,p):
    dimers = count_dimer(fasta_file,p)
    amino_acids = 'ACDEFGHIKLMNPQRSTVWXY'
    dimers_fq = dict()

    # this is to remove dimers which contains letters other than these 20 amino_acids,
    # like *
    for a1 in amino_acids:
        for a2 in amino_acids:
            dimers_fq[a1+a2] = dimers.get(a1+a2,0.0)
    number_of_aa_in_fasta = sum(dimers_fq.values())
    for key,value in dimers_fq.items(): dimers_fq[key] = value/number_of_aa_in_fasta
    return dimers_fq

def predict(fasta_file,model,means,stds,features,p):
    dimers_fq = get_dimer_frequency(fasta_file,p)

    Xs = list()
    for fea in features:
        Xs.append((dimers_fq[fea]-means[fea])/stds[fea])

    Xs = np.array(Xs).reshape([1,len(Xs)])

    pred_ogt = model.predict(Xs)[0]
    return np.around(pred_ogt,decimals=2)


def predOGT(args):
    infile = args.get('fasta', None)
    indir = args.get('indir', None)

    if args.get('o', None) is None: outf = sys.stdout
    else: outf = open(args['o'], 'w', encoding='utf-8')

    p = int(args.get('p',1))

    model, means, stds, features = load_model()
    outf.write('FileName\tpredOGT (C)\n')

    if infile is not None:
        pred_ogt = predict(infile,model,means,stds,features,p)
        # 使用 os.path.basename 替代 split('/') 以支持跨平台路径
        filename = os.path.basename(infile)
        outf.write('{0}\t{1}\n'.format(filename, pred_ogt))

    else:
        for name in os.listdir(indir):
            if name.startswith('.'): continue
            if not name.endswith('.fasta'): continue
            pred_ogt = predict(os.path.join(indir,name),model,means,stds,features,p)
            outf.write('{0}\t{1}\n'.format(name, pred_ogt))
    
    # 确保文件正确关闭
    if args.get('o', None) is not None:
        outf.close()

################################################################################




# Find homologues for a given enzyme with the same ec number
################################################################################

# to be completed
# external_data/

def download_external_data(link):
    """跨平台下载外部数据"""
    realpath = os.path.dirname(os.path.realpath(__file__))
    external_data_path = os.path.join(realpath, 'external_data')
    
    if not os.path.exists(external_data_path):
        os.makedirs(external_data_path)
    
    file_name = link.split('/')[-1]
    file_path = os.path.join(external_data_path, file_name)
    
    print_out('Downloading data from {0}'.format(link))
    
    try:
        # 使用 requests 进行跨平台下载
        response = requests.get(link, stream=True)
        response.raise_for_status()
        
        with open(file_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        
        print_out('Successfully downloaded {0}'.format(file_name))
    except Exception as e:
        print_out('Error downloading {0}: {1}'.format(file_name, str(e)))
        # 如果 requests 失败，尝试使用系统命令
        try:
            if platform.system() == 'Windows':
                # Windows 使用 PowerShell
                cmd = ['powershell', '-Command', 
                       f'Invoke-WebRequest -Uri "{link}" -OutFile "{file_path}"']
            else:
                # Unix-like 系统尝试 curl 或 wget
                if shutil.which('curl'):
                    cmd = ['curl', link, '-o', file_path]
                elif shutil.which('wget'):
                    cmd = ['wget', link, '-P', external_data_path]
                else:
                    raise Exception('No download tool available')
            
            subprocess.call(cmd)
        except Exception as e2:
            print_out('Failed to download using system commands: {0}'.format(str(e2)))
            raise


def get_enzymes_of_ec(ec,dfanno,temps,data_type):
    df = dfanno
    subdf = df.loc[ec,:]
    print_out('{0} sequences were found for {1}'.format(subdf.shape[0],ec))

    for i in range(len(df.columns)):
        if data_type.lower() == df.columns[i]: col_ind = i

    data = subdf.values
    data = data[data[:,col_ind]>temps[0],:]
    data = data[data[:,col_ind]<temps[1],:]

    index = pd.Index(data[:,0],name=df.columns[0])
    subdf = pd.DataFrame(data=data[:,1:],columns=df.columns[1:],index=index)
    print_out('{0} sequences were found between {1}'.format(subdf.shape[0],temps))

    return subdf

def build_fasta_for_given_ec(ec,uniprot_ids,brenda_seq_file,outdir):
    is_target = dict()
    for id in uniprot_ids: is_target[id] = True

    outfafile = os.path.join(outdir,'{0}_all.fasta'.format(ec))
    fhand = open(outfafile,'w')

    subseqs,ids = [],[]
    for rec in SeqIO.parse(brenda_seq_file,'fasta'):
        if not is_target.get(rec.id,False): continue
        fhand.write('>{0}\n{1}\n'.format(rec.id,rec.seq))
        subseqs.append(rec.seq)
        ids.append(rec.id)

    fhand.close()
    dfseqs = pd.DataFrame(data={'sequence':subseqs},index=ids)
    return dfseqs

def run_blastp(ec, seqfile, cpu_num, outdir, evalue):
    """运行 BLASTP，如果没有安装则跳过"""
    dbseq = os.path.join(outdir, '{0}_all.fasta'.format(ec))
    db = os.path.join(outdir, 'db')
    out = os.path.join(outdir, 'blast_{}.tsv'.format(ec))
    
    # 检查是否安装了 BLAST+
    if not shutil.which('makeblastdb') or not shutil.which('blastp'):
        print_out('Warning: BLAST+ not found. Please install NCBI BLAST+ to use sequence comparison features.')
        print_out('You can download it from: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download')
        return False
    
    try:
        # 创建数据库
        subprocess.run(['makeblastdb', '-dbtype', 'prot', '-in', dbseq, '-out', db], 
                      check=True, capture_output=True)
        
        # 运行 BLASTP
        subprocess.run(['blastp', '-query', seqfile, '-db', db, '-outfmt', '6', 
                       '-num_threads', str(cpu_num), '-out', out, '-evalue', evalue, 
                       '-max_hsps', '1'], check=True, capture_output=True)
        
        return True
    except subprocess.CalledProcessError as e:
        print_out('Error running BLAST: {0}'.format(str(e)))
        return False

def parse_blastp_results(outdir, ec):
    blastRes = dict()
    blastfile = os.path.join(outdir, 'blast_{}.tsv'.format(ec))
    fastafile = os.path.join(outdir, '{}_all.fasta'.format(ec))
    
    if not os.path.exists(blastfile):
        print_out('BLAST results file not found. Skipping sequence comparison.')
        return blastRes
    
    seqs = SeqIO.to_dict(SeqIO.parse(fastafile, 'fasta'))
    
    for line in open(blastfile):
        cont = line.split()
        target = cont[1]
        ident = float(cont[2])
        seq = seqs[target].seq
        cov = float(cont[3])/len(seq)*100
        
        blastRes[target] = (ident, cov, seq)
    
    # 跨平台删除临时文件
    try:
        if os.path.exists(blastfile):
            os.remove(blastfile)
        if os.path.exists(fastafile):
            os.remove(fastafile)
    except Exception as e:
        print_out('Warning: Could not remove temporary files: {0}'.format(str(e)))
    
    return blastRes

def get_info_for_selected_seqs(dfanno,uniprot_ids,ec):
    df = dfanno
    subdf = df.loc[ec,:]
    seqInfo = dict()
    # seqInfo = {uniprot_id:[domain,..,topt_source]}
    data = subdf.values
    for i in range(data.shape[0]):
        id = data[i,0]
        seqInfo[id] = [data[i,j] for j in range(data.shape[1])]
    return seqInfo


def build_output(blastRes,seqInfo,outdir,seqfile,dfanno):
    # two ouput files
    # 1. a fasta file containing all target sequence plus query
    # 2. a excel file containts the information of the target

    query = SeqIO.to_dict(SeqIO.parse(seqfile,'fasta'))
    query_id = list(query.keys())[0]
    query_seq = query[query_id].seq

    # write the fasta file
    outfasta = os.path.join(outdir,query_id+'_homologs.fasta')
    fhand = open(outfasta,'w')
    fhand.write('>{0}\n{1}\n'.format(query_id,query_seq))
    for id, rec in blastRes.items(): fhand.write('>{0}\n{1}\n'.format(id,rec[-1]))

    # build a dataframe and export it to excel file
    outcsv = os.path.join(outdir,query_id+'_homologs.tsv')
    #outexcel = os.path.join(outdir,query_id+'_homologs.xlsx')

    data = dict()
    cols = list(dfanno.columns) + ['identity(%)','coverage(%)','sequence']
    #cols = ['id','identity(%)','coverage(%)','domain','organism','source','growth_temp','sequence']
    # first line is for query
    data['uniprot_id'] = ['query']
    for col in cols[1:-1]: data[col] = [None]
    data['sequence'] = [query_seq]

    for id,rec in blastRes.items():
        for i in range(len(dfanno.columns)): data[dfanno.columns[i]] += [seqInfo[id][i]]

        data['identity(%)'] += [rec[0]]
        data['coverage(%)'] += [rec[1]]
        data['sequence'] += [rec[2]]

    df = pd.DataFrame(data=data,columns=cols)
    df.to_csv(outcsv,sep='\t')

    #writer = pd.ExcelWriter(outexcel)
    #df.to_excel(writer,'Sheet1')
    #writer.save()


def getEnzymes(args):

    fasta_link = 'https://zenodo.org/record/2539114/files/brenda_sequences_20180109.fasta'
    anno_link = 'https://zenodo.org/record/2539114/files/enzyme_ogt_topt.tsv'

    path = os.path.dirname(os.path.realpath(__file__))
    ext_dir = os.path.join(path, 'external_data')
    if not os.path.exists(ext_dir): 
        os.makedirs(ext_dir)
    
    seqfile = args.get('seq', None)
    ec = args.get('ec', None)
    outdir = args.get('outdir', None)
    cpu_num = int(args.get('p', 1))
    data_type = args.get('data_type', 'Topt')
    
    if data_type not in ['OGT', 'Topt']:
        sys.exit('Please check your --data_type option. Only OGT and Topt are supported')
    
    try: 
        temps = args['temp_range']
    except:
        sys.exit('Please specify --temp_range')
    
    temps = [float(item) for item in temps.split(',')]
    evalue = args.get('evalue', '1e-10')
    
    if ec is None: 
        sys.exit('Error: Please specify ec number.')
    
    if outdir is None: 
        outdir = './'
    if not os.path.exists(outdir): 
        os.makedirs(outdir)
    
    annofile = os.path.join(ext_dir, 'enzyme_ogt_topt.tsv')
    brenda_seq_file = os.path.join(ext_dir, 'brenda_sequences_20180109.fasta')
    
    if not os.path.isfile(annofile): 
        download_external_data(anno_link)
    if not os.path.isfile(brenda_seq_file): 
        download_external_data(fasta_link)
    
    dfanno = pd.read_csv(annofile, index_col=0, sep='\t')
    print_out('step 1: get all uniprot ids with the given ec number')
    subdf = get_enzymes_of_ec(ec, dfanno, temps, data_type)
    uniprot_ids = list(subdf.index)
    
    print_out('')
    print_out('step 2: get sequences and saving sequences to fasta format')
    dfseqs = build_fasta_for_given_ec(ec, uniprot_ids, brenda_seq_file, outdir)
    
    dfec_out = pd.merge(subdf, dfseqs, left_index=True, right_index=True, how='inner')
    
    if seqfile is None:
        subout = os.path.join(outdir, '{}_all.tsv'.format(ec))
        print_out('Saving results to tab-seperated format')
        dfec_out.to_csv(subout, sep='\t')
        sys.exit('Done!')
    
    print_out('step 3: run blastp')
    blast_success = run_blastp(ec, seqfile, cpu_num, outdir, evalue)
    
    if blast_success:
        print_out('step 4: get info of hits')
        seqInfo = get_info_for_selected_seqs(dfanno, uniprot_ids, ec)
        blastRes = parse_blastp_results(outdir, ec)
        print_out('{0} homologues were found by blast'.format(len(blastRes)))
        
        print_out('step 5: save results')
        build_output(blastRes, seqInfo, outdir, seqfile, dfanno)
    else:
        print_out('BLAST comparison skipped. Only basic enzyme information will be saved.')
        subout = os.path.join(outdir, '{}_all.tsv'.format(ec))
        dfec_out.to_csv(subout, sep='\t')
    
    print_out('Done!')
    
    # 清理数据库文件
    dbfile = os.path.join(outdir, 'db')
    try:
        for file in os.listdir(outdir):
            if file.startswith('db.'):
                os.remove(os.path.join(outdir, file))
    except Exception as e:
        print_out('Warning: Could not remove database files: {0}'.format(str(e)))

def main():
    args = parse_args()

    if args.get('help',None) is not None and args.get('method',None) is None:
        help_msg = '''
        Tome (Temperature optima for microorganisms and enzymes) is an open
        source suite for two purposes:
        (1) predict the optimal growth temperature from proteome sequences
        (2) get homologue enzymes for a given ec number with/without a sequence

        Tome Version 1.1.0 (built on Dec 2024)

        Main tools:
            predOGT     Predict optimal growth temperature(s) for one/many microorganisms
            getEnzymes  Get homologue enzymes for a given ec number with/without a sequence

        A detailed list of options can be obtained by calling 'tome predOGT --help'for
        predOGT or 'tome getEnzymes --help' for getEnzymes

        Gang Li
        2018-11-23
        '''
        sys.exit(help_msg)

    if args['method'] == 'predOGT':
        # check arguments and run the method
        all_args = ['fasta','indir','o','p','method','h','help','train']
        for arg in args.keys():
            if arg not in all_args:
                sys.exit('Unknown argument: {0}'.format(arg))

        if args.get('help',None) is not None or len(args)<1:
            help_msg = '''
        Usage:
        tome predOGT [Options]

        Options:
            --fasta input fasta file containing all protein sequence of a proteome.
                    Required for the prediction of OGT for one microorganism
            --indir directory that contains a list of fasta files. Each fasta file
                    is a proteome. Required for the prediction of OGT for a list of
                    microorganisms. Important: Fasta file names much endswith .fasta
            --train train the model again.
            -o      out file name, default: print on the terminal
            -p      number of threads, default is 1. if set to 0, it will use all cpus
                    available.

            '''
            sys.exit(help_msg)
        elif args.get('train',None) is not None: train_model()
        else: predOGT(args)

    elif args['method'] == 'getEnzymes':
        all_args = ['ec','seq','temp_range','outdir','p','method','h','help','evalue','data_type']
        for arg in args.keys():
            if arg not in all_args:
                sys.exit('Unknown argument: {0}'.format(arg))

        if args.get('-help',None) is not None:
            help_msg = '''
        Usage:
        tome getEnzymes [Options]

        Options:
            --ec         EC number. Required. 1.1.1.1, for instance
            --temp_range the temperature range that target enzymes should be in.
                         For example: 50,100. 50 is lower bound and 100 is upper
                         bound of the temperature.
            --data_type  [OGT or Topt], If OGT, Tome will find enzymes whose OGT
                         of its source organims fall into the temperature range
                         specified by --temp_range. If Topt, Tome will find enzymes
                         whose Topt  fall into the temperature range specified
                         by --temp_range. Default is Topt
            --seq        input fasta file which contains the sequence of the query
                         enzyme. Optional
            --outdir     directory for ouput files. Default is current working folder.
            -p           number of threads, default is 1. if set to 0, it will use all cpus
                         available.
            --evalue     evalue used in ncbi blastp. Default is 1e-10

            '''
            sys.exit(help_msg)

        else:getEnzymes(args)

    else:
        sys.exit('''Error:
        Please check your inputs a gain. Tome currently only provides preOGT and
        getEnzymes.
        \n''')

if __name__ == "__main__":
    main()
