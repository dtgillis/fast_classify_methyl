

import os
import numpy as np
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
from scipy.stats import pearsonr

data = rpy2.robjects.r['data']
smac = importr('smac')

dt = data('ex1.data')

data_ex1 = rpy2.robjects.r['ex1.data']
x = np.array(data_ex1.rx('ex1.x')).reshape((100, 2))
y = np.array(data_ex1.rx('ex1.y'))
x_new = np.array(data_ex1.rx('ex1.new.x')).reshape((100, 2))
y_new = np.array(data_ex1.rx('ex1.new.y')).reshape((100,))
fit = smac.cvsmac(x, y.T)

best_lambdas =  fit.rx('best.lambda')[0]
for best_lambda in best_lambdas:
    print float(best_lambda)
predictor = smac.predict_cvsmac(fit, x_new)
best_y = predictor.rx('best.pred.y')[0][0]
min_error = fit.rx('min.error')
print min_error[0][0]
#print fit.rx('model')









