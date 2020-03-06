import numpy as np

eps = 2.2204e-16

def filter_genes_zero(data):
	n_samples, n_genes = data.shape[0], data.shape[1]
	gene_nonzero = [False if (np.unique(data[:,col])==[0]).all() else True for col in range(n_genes)]
	return data[:, gene_nonzero]

def matrixNormalize(x):
	# normalize cells
	cells_norm2 = np.linalg.norm(x, axis=0)
	return x/cells_norm2, cells_norm2

def l_gene_select(ssc_score, pear_score, spear_score, cos_score):
	score_set = [ssc_score, pear_score, spear_score, cos_score]
	gene_inter, gene_inter_num = [], []
	for i in range(4):
		# ascend sort 
		score, sort_ind = np.sort(score_set[i]), np.argsort(score_set[i])
		# descend sort
		score, sort_ind = score[::-1], sort_ind[::-1]
		gene_num = len(score_set[i])
		thresh1 = int(np.round(0.1 * gene_num))
		thresh2 = int(np.round(0.5 * gene_num))

		gene_var = np.zeros((thresh2+1,))
		for j in np.arange(thresh1, thresh2+1):
			score1 = score[:j]
			score2 = score[j:]
			var1 = score1.var()
			var2 = score2.var()
			gene_var[j] = var1 + var2
		gene_var[:thresh1] = np.inf
		select_index = np.argmin(gene_var)
		gene_inter.append(sort_ind[:select_index])
		gene_inter_num.append(select_index)

	gene_select = list(set(gene_inter[0]) & set(gene_inter[1]) & set(gene_inter[2]) & set(gene_inter[3]))
	return gene_select, gene_inter, gene_inter_num

def l_gene_select_combine(ssc_score, pear_score, spear_score, cos_score):
	score_set = [ssc_score, pear_score, spear_score, cos_score]
	gene_inter, gene_inter_num = [], []
	for i in range(4):
		# ascend sort 
		score, sort_ind = np.sort(score_set[i]), np.argsort(score_set[i])
		# descend sort
		score, sort_ind = score[::-1], sort_ind[::-1]
		gene_num = len(score_set[i])
		thresh1 = int(np.round(0.1 * gene_num))
		thresh2 = int(np.round(0.5 * gene_num))

		gene_var = np.zeros((thresh2+1,))
		for j in np.arange(thresh1, thresh2+1):
			score1 = score[:j]
			score2 = score[j:]
			var1 = score1.var()
			var2 = score2.var()
			gene_var[j] = var1 + var2
		gene_var[:thresh1] = np.inf
		select_index = np.argmin(gene_var)
		gene_inter.append(sort_ind[:select_index])
		gene_inter_num.append(select_index)

	gene_slect_ssc_pear = np.intersect1d(gene_inter[0], gene_inter[1]) 
	gene_slect_ssc_spear = np.intersect1d(gene_inter[0], gene_inter[2])        
	gene_slect_ssc_cos = np.intersect1d(gene_inter[0], gene_inter[3]) 
	# combine
	gene_slect_combine = {}
	gene_slect_combine[0] = gene_slect_ssc_pear
	gene_slect_combine[1] = gene_slect_ssc_spear
	gene_slect_combine[2] = gene_slect_ssc_cos
	# combine 
	gene_slect_combine[3] = np.intersect1d(gene_slect_ssc_pear, gene_inter[2]) 
	gene_slect_combine[4] = np.intersect1d(gene_slect_ssc_pear, gene_inter[3])     
	gene_slect_combine[5] = np.intersect1d(gene_slect_ssc_spear, gene_inter[3])    

	gene_select = np.intersect1d(gene_slect_combine[3], gene_inter[3])
	return gene_select, gene_slect_combine, gene_inter_num

def l_choose_edge(pear_sim, spear_sim, cos_sim):
	m, n = pear_sim.shape[0], pear_sim.shape[1]
	print(m, n)
	if m > 5000:
		edge_num = 100
	else:
		edge_num = int(np.round(m*0.1))

	index1 = np.argsort(pear_sim, axis=1)
	pear_index = index1[:, -edge_num:]

	index2 = np.argsort(spear_sim, axis=1)
	spear_index = index2[:, -edge_num:]

	index3 = np.argsort(cos_sim, axis=1)
	cos_index = index3[:, -edge_num:]

	elected_edge = []
	for i in range(m):
		tmp = list(set(pear_index[i]) | set(spear_index[i]) | set(cos_index[i]))
		tmp = np.sort(tmp)
		elected_edge.append(tmp)
	return elected_edge

def l_choose_edge_combine(pear_sim, spear_sim, cos_sim):
	m, n = pear_sim.shape[0], pear_sim.shape[1]
	# print(m, n)
	if m > 5000:
		edge_num = 100
	else:
		edge_num = int(np.round(m*0.1))

	index1 = np.argsort(pear_sim, axis=1)
	pear_index = index1[:, -edge_num:]

	index2 = np.argsort(spear_sim, axis=1)
	spear_index = index2[:, -edge_num:]

	index3 = np.argsort(cos_sim, axis=1)
	cos_index = index3[:, -edge_num:]
	# selected edge with different combination of similarities
	selected_edge_combine = {}
	selected_edge_combine[0] = pear_index
	selected_edge_combine[1] = spear_index
	selected_edge_combine[2] = cos_index

	elected_edge = []
	for i in range(m):
		tmp = np.union1d(pear_index[i,:], spear_index[i,:])
		elected_edge.append(np.sort(tmp))
	selected_edge_combine[3] = elected_edge

	elected_edge = []
	for i in range(m):
		tmp = np.union1d(pear_index[i,:], cos_index[i,:])
		elected_edge.append(np.sort(tmp))
	selected_edge_combine[4] = elected_edge

	elected_edge = []
	for i in range(m):
		tmp = np.union1d(spear_index[i,:], cos_index[i,:])
		elected_edge.append(np.sort(tmp))
	selected_edge_combine[5] = elected_edge
	return selected_edge_combine

def l_enhance(data, select_edge):
	data = np.abs(data)
	m, n = data.shape[0], data.shape[1]
	test_data = data.copy()

	RA_score = np.zeros((m,n))
	WRA_score = np.zeros((m,n))

	for i in range(m):
		edge_len = len(select_edge[i])
		for j in range(edge_len):
			if data[i, select_edge[i][j]] == 0:
				RA, WRA = 0, 0
				for z in range(m):
					if data[i,z]!=0 and data[select_edge[i][j],z]!=0:
						neighbor_num = list(data[z]!=0).count(True)
						if neighbor_num==0:
							print('a error happened')
						else:
							RA = RA + 1/neighbor_num
							WRA = WRA + (data[i,z] + data[select_edge[i][j], z])/neighbor_num
				RA_score[i, select_edge[i][j]] = RA
				WRA_score[i, select_edge[i][j]] = WRA


	WRA_value = WRA_score + WRA_score.T
	test_data0 = data + WRA_value
	test_data = test_data0 - np.diag(test_data0.diagonal())
	return test_data

def Estimate_Number_of_Clusters_Gap(W, Num):
	N = W.shape[0]
	DN = np.diag(1./np.sqrt(W.sum(axis=0)+eps))
	LapN = np.diag(np.ones(N,)) - DN.dot(W).dot(DN)

	_, U, _ =  np.linalg.svd(LapN)
	eigvalue1 = U.diagonal()
	eig1_select = eigvalue1[-Num:]
	eigvalue2 = U.diagonal()
	eig2_select = eigvalue2[-(Num+1):-1]

	tmp = np.abs(eig2_select - eig1_select)
	# ascend then reverse
	ind = np.argsort(tmp)[::-1]
	K1 = Num - ind[0] + 1
	return K1

