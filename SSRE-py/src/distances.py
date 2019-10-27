from sklearn.metrics.pairwise import pairwise_distances
import numpy as np
import pandas as pd

def cosine(data):
	cosine = 1 - pairwise_distances(data, metric='cosine')
	return np.where((cosine<=1) & (cosine>=0), cosine, np.zeros(cosine.shape))

def euclidean(data):
	dist = 1 - pairwise_distances(data, metric='euclidean')
	# dist transform sim
	# sim = np.exp(-dist/dist.max())
	return np.where((dist<=1)&(dist>=0), dist, np.zeros(dist.shape))

def pearson(data):
	df = pd.DataFrame(data.T)
	pear_ = df.corr(method='pearson')
	return np.where(pear_>=0, pear_, np.zeros(shape=(pear_.shape)))

def spearman(data):
	df = pd.DataFrame(data.T)
	spear_ = df.corr(method='spearman')
	return np.where(spear_>=0, spear_, np.zeros(shape=(spear_.shape)))

def all_similarities(data):
	return pearson(data), spearman(data), cosine(data), euclidean(data)
