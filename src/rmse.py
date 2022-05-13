import numpy as np

def rmse(calc, true):
    '''
    Calculate the root mean square error (RMSE)

    Args:
        calc (array): simulate value
        true (array): true value

    Returns:
        rmse (float): root mean square error
    '''
    assert calc.shape == true.shape, "Dimension of comparison should be the same, %s, %s" % (calc.shape, true.shape)
    if len(true) > 1.:
        rmse = np.zeros(len(true) - 1)
        for i in range(1, len(true)):
            rmse[i-1] = np.sqrt(np.average(np.square(calc[i] - true[i])))
    else:
        rmse = np.sqrt(np.mean(np.square(calc - true)))
    return rmse
