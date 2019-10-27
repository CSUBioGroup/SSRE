import numpy as np
import matplotlib.pyplot as plt

realmin = 2.2251e-308
def tsne_l(P, label=None, n_dims=2):
	#P is a similarity matrix
	# Check initial solution
	if not isinstance(n_dims, int) and len(np.reshape(n_dims, -1))>1:
		initial_solution = True
		ydata = n_dims
		n_dims = y_data.shape[1]
	else:
		initial_solution = False
	P = P / np.sum(P)
	n = P.shape[0]
	momentum = .08;                                     # initial momentum
	final_momentum = .1;                               # value to which momentum is changed
	mom_switch_iter = 250;                              # iteration at which momentum is changed
	stop_lying_iter = 100;                              # iteration at which lying about P-values is stopped
	max_iter = 1000;                                    # maximum number of iterations
	epsilon = 500;                                      # initial learning rate
	min_gain = .01;   

	for i in range(n):
		P[i,i] = 0.
	P = 0.5*(P+P.T)
	P /= P.sum()
	P = np.where(P>=realmin, P, np.ones(P.shape)*realmin)
	const = (P*np.log(P)).sum()
	if not initial_solution:
		P = P * 4

	if not initial_solution:
		ydata = 0.0001 * np.random.randn(n, n_dims)

	y_incs = np.zeros(ydata.shape)
	gains = np.ones(ydata.shape)

	for itera in range(max_iter):
		sum_ydata = (ydata**2).sum(axis=1)

		# column add
		denomi = sum_ydata.T.reshape((-1,1)) - 2* ydata.dot(ydata.T)
		# row add
		denomi = denomi + sum_ydata
		denomi += 1
		num = 1/denomi 

		# reset the diagal
		for i in range(num.shape[0]):
			num[i,i] = 0
		Q = num / num.sum()
		Q = np.where(Q>=realmin, Q, np.ones(Q.shape)*realmin)

		L = (P-Q) * num
		y_grads = 4 * (np.diag(L.sum(axis=0)) - L).dot(ydata)

		gains = (gains + 0.2) * (np.sign(y_grads) != np.sign(y_incs)) + (gains * 0.8) * (np.sign(y_grads)==np.sign(y_incs))
		gains = np.where(gains>=min_gain, gains, np.ones(gains.shape)*min_gain)
		y_incs = momentum * y_incs - epsilon * (gains * y_grads)
		ydata = ydata + y_incs
		ydata = ydata - ydata.mean(axis=0)
		ydata = np.where(ydata>=-100, ydata, np.ones(ydata.shape)*(-100))
		ydata = np.where(ydata<=100, ydata, np.ones(ydata.shape)*100)

		if itera == mom_switch_iter:
			momentum = final_momentum
		if itera == stop_lying_iter and not initial_solution:
			P = P / 4

	return ydata
