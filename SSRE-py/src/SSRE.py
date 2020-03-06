from sklearn.metrics import normalized_mutual_info_score, adjusted_rand_score
from time import time
import numpy as np

import os
from os.path import join, abspath, dirname
import sys
sys.path.append(dirname(abspath(__file__)))

from laplacian import LaplacianScore, spectralclustering
from distances import all_similarities
from utils import filter_genes_zero, matrixNormalize, l_gene_select_combine, l_enhance, l_choose_edge_combine, l_gene_select,l_choose_edge
from admm_lasso import admmLasso_mat_func

eps = 2.2204e-16

def SSRE(data, label, omega):
	# calculate pairwise similarity
	st = time()
	n_clusters = len(set(label))
	data = filter_genes_zero(data)  # row is a cell
	pairwise_data0 = data
	pear_sim0,spear_sim0,cos_sim0, _ = all_similarities(pairwise_data0)
	ed = time()
	print('phase1, ', ed-st)

	# gene score
	st = time()
	pear_score  = LaplacianScore(pairwise_data0, pear_sim0);
	spear_score = LaplacianScore(pairwise_data0, spear_sim0);
	cos_score   = LaplacianScore(pairwise_data0, cos_sim0);
	ed = time()
	print('phase2, ', ed-st)

	# calculate sparse similarity
	st = time()
	ssc_data0, _ = matrixNormalize(data.T)	# (features, samples)
	CMat3 = admmLasso_mat_func(ssc_data0, False, omega)
	C0 = np.abs(CMat3) + np.abs(CMat3.T) + eps
	ssc_score = LaplacianScore(ssc_data0.T, C0)
	ed = time()
	print('phase3, ', ed-st)

	# gene selection
	st = time()
	gene_select, _ , _ = l_gene_select(ssc_score, pear_score, spear_score, cos_score)
	gene_select1 = np.sort(gene_select)
	data_select = data[:, gene_select1]
	ed = time()
	print('phase4, ', ed-st)

	# calculate sparse similarity
	st = time()
	ssc_data1, _ = matrixNormalize(data_select.T)
	CMat3 = admmLasso_mat_func(ssc_data1, False, omega)
	C1 = np.abs(CMat3) + np.abs(CMat3.T) + eps
	ed = time()
	print('phase5, ', ed-st)

	# calculate pairwise similarity
	st = time()
	pairwise_data = data_select
	[pear_sim1, spear_sim1, cos_sim1, _] = all_similarities(pairwise_data)

	elected_edges = l_choose_edge(pear_sim1, spear_sim1, cos_sim1)

	CC1 = C1 - eps
	test_data = l_enhance(CC1, elected_edges)
	ed = time()
	print('phase6, ', ed-st)

	st = time()
	clf = spectralclustering(test_data, n_clusters)

	NMI, ARI = normalized_mutual_info_score(label, clf), adjusted_rand_score(label, clf)
	ed = time()
	print('phase7, ', ed-st)

	print('Result, NMI: {:.4f}, ARI: {:.4f}'.format(NMI, ARI))
	return NMI, ARI, clf, test_data


def SSRE_combine(data, label, omega):
	# calculate pairwise similarity
	st = time()
	n_clusters = len(set(label))
	data = filter_genes_zero(data)  # row is a cell
	pairwise_data0 = data
	pear_sim0,spear_sim0,cos_sim0, _ = all_similarities(pairwise_data0)
	ed = time()
	print('step1, {}s'.format(ed-st))

	# gene score
	st = time()
	pear_score  = LaplacianScore(pairwise_data0, pear_sim0);
	spear_score = LaplacianScore(pairwise_data0, spear_sim0);
	cos_score   = LaplacianScore(pairwise_data0, cos_sim0);
	ed = time()
	print('step2, {}s'.format(ed-st))

	# calculate sparse similarity
	st = time()
	ssc_data0, _ = matrixNormalize(data.T)	# (features, samples)
	CMat3 = admmLasso_mat_func(ssc_data0, False, omega)
	C0 = np.abs(CMat3) + np.abs(CMat3.T)
	ssc_score = LaplacianScore(ssc_data0.T, C0)
	ed = time()
	print('step3, {}s'.format(ed-st))

	# gene selection with different combination
	gene_select, gene_slect_combine, _ = l_gene_select_combine(ssc_score, pear_score, spear_score, cos_score)
	
	NMI, ARI = np.zeros(6,), np.zeros(6,)
	for i in range(6):
		gene_select1 = np.sort(gene_slect_combine[i])	
		data_select = data[:, gene_select1]
		
		# calculate sparse similarity
		ssc_data1, _ = matrixNormalize(data_select.T)
		CMat3 = admmLasso_mat_func(ssc_data1, False, omega)
		C1 = np.abs(CMat3) + np.abs(CMat3.T) 
		# print('enhanced_data---------------------------------------------')
		# print(C1[:10,:16])
		# calculate pairwise similarity
		pairwise_data = data_select
		pear_sim1, spear_sim1, cos_sim1, _ = all_similarities(pairwise_data)

		selected_edge_combine = l_choose_edge_combine(pear_sim1, spear_sim1, cos_sim1)
		elected_edges = selected_edge_combine[i]
		# CC1 = C1 - eps
		# test_data = l_enhance(CC1, elected_edges)
		test_data = l_enhance(C1, elected_edges)
		cluster = spectralclustering(test_data, n_clusters)
		NMI[i] = normalized_mutual_info_score(label, cluster)
		ARI[i] = adjusted_rand_score(label, cluster)
		#print('Result, NMI: {:.4f}, ARI: {:.4f}'.format(NMI, ARI))
	return NMI, ARI
