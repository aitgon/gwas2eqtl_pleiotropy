import scipy
from matplotlib import pyplot as plt

from gwas2eqtl_pleiotropy.constants import label_fontsize


def boxenplot_with_mannwhitneyu(group1, group2, x1, x2, y, h):
    col = 'k'
    statistic, pvalue = scipy.stats.mannwhitneyu(group1, group2)
    plt.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c=col)
    if pvalue < 1.00e-04:
        pvalue_asterisk = "****"
    elif pvalue < 1.00e-03:
        pvalue_asterisk = "***"
    elif pvalue < 1.00e-02:
        pvalue_asterisk = "**"
    elif pvalue < 5.00e-02:
        pvalue_asterisk = "*"
    else:
        pvalue_asterisk = "ns"
    plt.text((x1 + x2) * .5, y + h, pvalue_asterisk, ha='center', va='bottom', color=col, fontsize=label_fontsize)
