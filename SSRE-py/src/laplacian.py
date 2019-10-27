import numpy as np
from sklearn.cluster import KMeans

eps = 2.2204e-16

def eps2C(x, c=10000):
	return np.where(x>=1e-12, x, c*np.ones(shape=x.shape))


def LaplacianScore(x, w):
	# x in (samples, features)
	n_samples, n_feat = x.shape[0], x.shape[1]

	if w.shape[0] != n_samples:
		raise Exception("W.shape not match X.shape")

	D = np.diag(np.sum(w, axis=1)) # (n_samples,)
	D2 = np.sum(w, axis=1) # (n_samples,)
	L = w

	tmp1 = (D2.T).dot(x)
	DPrime = np.sum((x.T.dot(D)).T * x, axis=0) - tmp1 * tmp1/np.sum(D2)
	LPrime = np.sum((x.T.dot(L)).T * x, axis=0) - tmp1 * tmp1/np.sum(D2)

	DPrime = eps2C(DPrime, c=10000)
	a1=np.sum(D)
	a2=np.sum((x.T.dot(D)).T * x, axis=0)
	a3=tmp1 * tmp1/np.sum(D)
	a4=(x.T.dot(D)).T * x
	a7=((x.T).dot(D)).T * x
	a5=tmp1 * tmp1
	a6=x.T.dot(D)
	a9=np.dot(x.T,D)

	Y = LPrime / DPrime
	#Y = Y.T#lzl edit
	return Y

def spectralclustering(data, n_clusters):
	N = data.shape[0]
	maxiter = 1000  # max iteration times
	replic = 100  # number of time kmeans will be run with diff centroids

	DN = np.diag(1/np.sqrt(np.sum(data, axis=0) + eps))
	lapN = np.eye(N) - DN.dot(data).dot(DN)
	U, A, V = np.linalg.svd(lapN)
	V = V.T
	kerN = V[:, N - n_clusters :N ]
	normN = np.sum(kerN**2, 1)**(0.5)
	kerNS = (kerN.T/(normN + eps)).T
	# kmeans
	clf = KMeans(n_clusters=n_clusters, max_iter=maxiter, n_init=replic)
	return clf.fit_predict(kerNS)