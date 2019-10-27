import scipy.io as sio 

def load_data(path):
	data = sio.loadmat(path)
	X = data['in_X']
	labels = data['true_labs'].reshape(-1)
	return X, labels