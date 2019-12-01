# example: https://stephens999.github.io/fiveMinuteStats/likelihood_ratio_simple_models.html
import numpy as np

data = np.array([1, 0, 1, 0, 0, 1])
f_S  = np.array([0.40, 0.12, 0.21, 0.12, 0.02, 0.32])
f_F  = np.array([0.80, 0.20, 0.11, 0.17, 0.23, 0.25])

# p (x|M) = ∏( f^x (1-f)^(1-x) )
def prob(f, x):
    foo = f**x * (1 - f)**(1 - x)
    return np.prod(foo)

likelihood_S = prob(f_S, data)
likelihood_F = prob(f_F, data)
likelihood_ratio = likelihood_S / likelihood_F
print('      LR: {}'.format(likelihood_ratio))
print('  Log LR: {}'.format(np.log(likelihood_ratio)))
print('Log10 LR: {}'.format(np.log10(likelihood_ratio)))
