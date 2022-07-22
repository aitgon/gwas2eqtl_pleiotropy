import itertools

import numpy
import pandas
from scipy.sparse import csr_matrix

#%%

count_per_rsid_gwas_tsv_path = "/home/gonzalez/Repositories/eqtl2gwas_pleiotropy/out/gwas413/genome/5e-08/1000000/cmpt_pleiotropic_regions.py/region_window_100000.tsv"
count_per_rsid_gwas_df = pandas.read_csv(count_per_rsid_gwas_tsv_path, sep="\t", header=0)

#%%
gwas_category_lst = count_per_rsid_gwas_df['gwas_category_lst']
gwas_category_2dlst = gwas_category_lst.apply(lambda x: x.split(','))

#%%
documents = gwas_category_2dlst
allowed_words = sorted(set([x for xs in documents for x in xs]))

#%%
def create_co_occurences_matrix(allowed_words, documents):
    print(f"allowed_words:\n{allowed_words}")
    print(f"documents:\n{documents}")
    word_to_id = dict(zip(allowed_words, range(len(allowed_words))))
    documents_as_ids = [numpy.sort([word_to_id[w] for w in doc if w in word_to_id]).astype('uint32') for doc in documents]
    row_ind, col_ind = zip(*itertools.chain(*[[(i, w) for w in doc] for i, doc in enumerate(documents_as_ids)]))
    data = numpy.ones(len(row_ind), dtype='uint32')  # use unsigned int for better memory utilization
    max_word_id = max(itertools.chain(*documents_as_ids)) + 1
    docs_words_matrix = csr_matrix((data, (row_ind, col_ind)), shape=(len(documents_as_ids), max_word_id))  # efficient arithmetic operations with CSR * CSR
    words_cooc_matrix = docs_words_matrix.T * docs_words_matrix  # multiplying docs_words_matrix with its transpose matrix would generate the co-occurences matrix
    words_cooc_matrix.setdiag(0)
    print(f"words_cooc_matrix:\n{words_cooc_matrix.todense()}")
    return words_cooc_matrix, word_to_id

#%%
words_cooc_matrix, word_to_id = create_co_occurences_matrix(allowed_words, documents)

#%%
countij = pandas.DataFrame.sparse.from_spmatrix(words_cooc_matrix, index=word_to_id.keys(), columns=word_to_id.keys())
countij = countij.loc[countij.index.sort_values(), ]
countij = countij[countij.columns.sort_values()]

#%% remove bottom triangular
countij_triu2 = pandas.DataFrame(numpy.triu(countij.values), index=countij.index, columns=countij.columns)

#%% total nb of pairs
countij_triu2_n = countij_triu2.sum().sum()

#%%
freqij2 = countij_triu2/countij_triu2_n

#%%
counti = countij.sum()
freqi = counti/(countij_triu2_n*2)

#%%
fc0 = freqij2.divide(freqi, axis=0)
fc1 = fc0.divide(freqi, axis=1)

#%%
numpy.log2(fc1).min()
