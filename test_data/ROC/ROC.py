import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import matthews_corrcoef, roc_curve, auc


def ACC(y_true, y_pred):
    x = 0
    for i in range(y_pred.shape[0]):
        a = round(y_true[i,])
        b = round(y_pred[i,])
        if (b == a):
            x = x + 1
    return round(x / y_pred.shape[0], 2)


def MCC(y_true, y_pred):
    x = matthews_corrcoef(y_true, y_pred)
    return x


def plotauc(drug, fpr, tpr, auc, odir):
    plt.figure()
    plt.title(drug)
    plt.plot(fpr, tpr, 'b', label="AUC = %0.2f" % auc)
    plt.legend(loc="lower right")
    plt.ylabel("True Positive Rate")
    plt.xlabel("False Positive Rate")
    plt.gcf().savefig(odir + "/roc." + str(auc) + ".pdf")
    return


data = pd.read_csv("ROC.csv", sep=",", index_col=False, header=0)


real = np.array(data['real'])
pre = np.array(data['predict'])
a = ACC(real, pre)
m = MCC(real, pre)
print("cal acc and mcc:{}, {}".format(a, m))
fpr, tpr, threshold = roc_curve(real, pre)
rocauc = auc(fpr, tpr)
od = "."
print("cal auc:{}".format(rocauc))
plotauc("ROC", fpr, tpr, rocauc, od)
