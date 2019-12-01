# https://stephens999.github.io/fiveMinuteStats/likelihood_function.html
import numpy as np
import matplotlib.pyplot as plt

q = np.linspace(0,1, 100)

def likelihood_fn(x):
    return np.power(x, 30) * np.power((1-x), 70)

l_val = likelihood_fn(q)
# print(q)
# print(l_val)
max_index = np.argmax(l_val)

textstr = '\n'.join(('L_max: {0:.3e}'.format(l_val[max_index]),
                     'q: {0:.4f}'.format(q[max_index])))
print(textstr)

if True:
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.plot(q, l_val)
    plt.axvline(x=q[max_index], color='r', ls='--')
    plt.text(q[max_index]+0.05, l_val[max_index], textstr, fontsize=12, verticalalignment='top', bbox=props)
    plt.xlabel('q (probability)')
    plt.ylabel('likelihood')
    plt.show()


# normalization
l_val_normalized = l_val / l_val[max_index]
textstr += '\nnorm: {} (of course!)'.format(l_val_normalized[max_index])

if True:
    plt.plot(q, l_val_normalized)
    plt.axvline(x=q[max_index], color='r', ls='--')
    plt.text(q[max_index]+0.05, l_val_normalized[max_index], textstr, fontsize=12, verticalalignment='top', bbox=props)
    plt.xlabel('q (probability)')
    plt.ylabel('L(q)/L(q_hat)')
    plt.show()

def log_likelihood_fn(x):
    return 30*np.log(x) + 70 * np.log(1-x)

log_l_val = log_likelihood_fn(q)
print(log_l_val)
# print(np.argmax(log_l_val), log_l_val[np.argmax(log_l_val)])
log_l_val_normalized = log_l_val - log_likelihood_fn(0.3)
if True:
    plt.plot(q, log_l_val_normalized)
    plt.xlabel('q')
    plt.ylabel('log L(q) / log L(q_hat)')
    plt.ylim(-10,0)
    plt.show()
