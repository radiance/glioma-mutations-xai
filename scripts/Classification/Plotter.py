import shap
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import random
import matplotlib.pyplot as pl

from sklearn.metrics import confusion_matrix


def plot_confusion_matrix(y_test, y_pred):
    cm = pd.DataFrame(confusion_matrix(y_test, y_pred), columns=["0-18", "19-70", "70+"], index=["0-18", "19-70", "70+"]) #columns=["0-23", "43-50", "50+"], index=["0-23", "43-50", "50+"]
    sns.heatmap(cm,annot=True,cmap='Blues', fmt='g')
    plt.show()
    #plt.savefig('heatmap.png')


def plot_feature_importance_for_class(model, X_train):
    explainer = shap.TreeExplainer(model, data=shap.sample(X_train, 100), feature_dependence="interventional")
    shap_values = explainer.shap_values(X_train)
    #shap.dependence_plot("rank(10)", shap_values, X_train) #raise TypeError("The passed shap_values are a list not an array! If you have a list of explanations try " \ # TypeError: The passed shap_values are a list not an array! If you have a list of explanations try passing shap_values[0] instead to explain the first output class of a multi-output model.
    shap.summary_plot(shap_values, X_train, plot_type="bar", class_names=model.classes_, color=pl.get_cmap("tab10")) #labels=model.classes_
    # IMPORTANT for some reason the three lines below might break if the line above isn't commmented out
    #shap.summary_plot(shap_values[0], X_train, class_names=model.classes_)
    #shap.summary_plot(shap_values[1], X_train, class_names=model.classes_)
    #shap.summary_plot(shap_values[2], X_train, class_names=model.classes_) #show=False
    #Feature values in pink cause to increase the prediction.
    #Size of the bar shows the magnitude of the feature's effect.
    #Feature values in blue cause to decrease the
    plt.show()


def plot_prediction_desicion(model, X_test, pred, row_idx):
    #The decision plot below shows the model’s multiple outputs for a single observation
    #the dashed line is the prediction of our classifier
    explainer = shap.TreeExplainer(model, data=shap.sample(X_test, 100), feature_dependence="interventional")
    shap_values = explainer.shap_values(X_test)
    shap.multioutput_decision_plot([1, 2, 3], shap_values,
                                   row_index=row_idx,
                                   feature_names=list(X_test.columns) ,
                                   highlight=int(pred[row_idx]),
                                   legend_labels=["0-18", "19-70", "70+"], #legend_labels=["0-23", "24-50", "50+"],
                                   legend_location='lower right')
    plt.show()

#def plot_data_balance(data_frame, label_col):
    #data_frame.groupby(label_col).plot.bar(ylim=0) #data_frame.groupby(label_col).Age.count().plot.bar(ylim=0)
    #plt.show()
