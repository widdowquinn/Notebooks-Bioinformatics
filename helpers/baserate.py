# baserate.py
#
# Helper code for the baserate_fallacy.ipynb notebook in the
# Notebooks-Bioinformatics repo

import matplotlib.pylab as pylab


def p_correct_given_pos(sens, fpr, b):
    """Returns a simple Bayesian probability for the probability
    that a prediction is correct, given that the prediction
    was positive, for the prevailing sensitivity (sens),
    false positive rate (fpr) and base rate of positive 
    examples.
    """
    assert 0 <= sens <= 1, "Sensitivity must be in range [0,1]"
    assert 0 <= fpr <= 1, "FPR must be in range [0,1]"
    return sens * b / (sens * b + fpr * (1 - b))


def plot_prob_effector(sens, fpr, xmax=1, baserate=0.1):
    """Plots a line graph of P(effector|positive test) against
    the baserate of effectors in the input set to the classifier.
        
    The baserate argument draws an annotation arrow
    indicating P(pos|+ve) at that baserate
    """
    assert 0.1 <= xmax <= 1, "Max x axis value must be in range [0,1]"
    assert 0.01 <= baserate <= 1, "Baserate annotation must be in range [0,1]"
    baserates = pylab.arange(0, 1.05, xmax * 0.005)  
    probs = [p_correct_given_pos(sens, fpr, b) for b in baserates]
    pylab.plot(baserates, probs, 'r')
    pylab.title("P(eff|pos) vs baserate; sens: %.2f, fpr: %.2f" % (sens, fpr))
    pylab.ylabel("P(effector|positive)")
    pylab.xlabel("effector baserate")
    pylab.xlim(0, xmax)
    pylab.ylim(0, 1)
    # Add annotation arrow
    xpos, ypos = (baserate, p_correct_given_pos(sens, fpr, baserate))
    if baserate < xmax:
        if xpos > 0.7 * xmax:
            xtextpos = 0.05 * xmax
        else:
            xtextpos = xpos + (xmax-xpos)/5.
        if ypos > 0.5:
            ytextpos = ypos - 0.05
        else:
            ytextpos = ypos + 0.05
        pylab.annotate('baserate: %.2f, P(pos|+ve): %.3f' % (xpos, ypos), 
                       xy=(xpos, ypos), 
                       xytext=(xtextpos, ytextpos),
                       arrowprops=dict(facecolor='black', shrink=0.05))
    else:
        pylab.text(0.05 * xmax, 0.95, 'baserate: %.2f, P(pos|+ve): %.3f' %
                   (xpos, ypos))
