import numpy as np
from os.path import join, abspath, dirname
import sys
sys.path.append(dirname(abspath(__file__)))

from preprocessing import matrixNormalize
from time import time

def computeLambda_mat(Y, P=None):
	if P is None:
		P = Y.copy()
	n_samples = Y.shape[1]
	T = P.T.dot(Y)
	T[:n_samples] = T[:n_samples] - np.diag(T[:n_samples].diagonal())
	lamda = np.min(np.max(T, 0))	
	return lamda

def errorCoef(Z, C):
	return np.abs(Z-C).max()

def errorLinSys(P, Z):
	R, N = Z.shape[0], Z.shape[1]
	if R>N:
		E = P[:, N+1:].dot(Z[N+1:])
		Y = P[:, :N]
		Y0 = Y - E
		C = Z[:N]
	else:
		Y = P
		Y0 = P
		C = Z
	Yn, norm = matrixNormalize(Y0)
	norm = norm.reshape((1, -1))
	M = np.tile(norm, (Y.shape[0], 1))
	S = Yn - Y.dot(C)/M
	err = np.sqrt(np.max(np.sum(S**2, axis=0)))
	return err

def admmLasso_mat_func(Y, affine=False, alpha=800, thr=2e-4, maxIter=200):
	if isinstance(alpha, int) or isinstance(thr, float):
		alpha = [alpha]
	if isinstance(thr, float) or isinstance(alpha, int):
		thr = [thr]
		
	alpha1, alpha2 = alpha[0], alpha[-1]
	thr1, thr2     = thr[0], thr[-1]
	n_samples = Y.shape[1]

	# setting penalty 
	mu1 = alpha1 * 1/computeLambda_mat(Y)
	mu2 = alpha2

	if not affine:
		# initialization
		A = np.linalg.inv(mu1*(Y.T.dot(Y)) + mu2*np.eye(n_samples))
		C1 = np.zeros(shape=(n_samples, n_samples))
		lamda2 = np.zeros(shape=(n_samples, n_samples))
		err1, err2 = [10 * thr1], [10 * thr2]
		i = 0
		start_stamp = time()
		while err1[i] > thr1 and i<maxIter:
			# updating z
			Z = A.dot(mu1*(Y.T.dot(Y)) + mu2*(C1-lamda2/mu2))
			Z = Z - np.diag(Z.diagonal()) 

			# updating c
			tmp = np.abs(Z+lamda2/mu2) - 1/mu2*np.ones(shape=(n_samples, n_samples))
			C2 = np.where(tmp>=0, tmp, np.zeros((tmp.shape))) * np.sign(Z+lamda2/mu2)
			C2 = C2 - np.diag(C2.diagonal())

			# updating lagrnge multipliers
			lamda2 += mu2 * (Z-C2)

			# computing errors
			err1.append(errorCoef(Z, C2))
			err2.append(errorLinSys(Y,Z))

			C1 = C2.copy() 
			i += 1
		end_stamp = time()
		print('lasso running time {}s, err1: {:.4f}, err2:{:.4f}, iter:{}'.format(end_stamp-start_stamp, err1[-1], err2[-1], i))
	else:
		# initialization
		A = np.linalg.inv(mu1*(Y.T.dot(Y)) + mu2*np.eye(n_samples) + mu2*np.ones((n_samples,n_samples)))
		C1 = np.zeros(shape=(n_samples, n_samples))
		lamda2 = np.zeros(shape=(n_samples, n_samples))
		lamda3 = np.zeros(shape=(1, n_samples))
		err1, err2, err3 = [10*thr1], [10*thr2], [10*thr1]
		i = 1
		while (err1[i]>thr1 or err3[i]>thr1) and i < maxIter:
			# updating z
			Z = A.dot(mu1*(Y.T.dot(Y)) + mu2*(C1-lamda2/mu2) + mu2*np.ones((n_samples,1)).dot(np.ones((1,n_samples))-lamda3/mu2))
			Z = Z - np.diag(Z.diagonal())

			# c
			tmp = np.abs(Z+lamda2/mu2) - 1/mu2*np.ones((n_samples,n_samples))
			C2 = np.max(tmp>=0, tmp, np.zeros(tmp.shape)) * np.sign(Z+lamda2/mu2)
			C2 = C2 - np.diag(C2.diagonal())

			# lagrange multipliers
			lamda2 += mu2*(Z-C2)
			lamda3 += mu2*(np.ones((1,n_samples)).dot(Z) - np.ones((1,n_samples)))

			# errors
			err1.append(errorCoef(Z, C2))
			err2.append(errorLinSys(Y, Z))
			err3.append(errorCoef(np.ones((1,n_samples)).dot(Z), np.ones((1, n_samples))))

			# 
			C1 = C2.copy()
			i += 1
		print('err1: {:.4f}, err2:{:.4f}, err3:{:.4f}, iter:{}'.format(err1[-1], err2[-1], err3[-1], i))
	return C2