import pandas as pd
import Plotter as Plotter

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC


class Classifier:
    data = pd.DataFrame()
    labels = pd.DataFrame()

    def __init__(self, labels, data):
        self.data = data
        self.labels = labels

    def predict(self):
        X_train, X_test, y_train, y_test = train_test_split(self.data, self.labels, test_size=0.20, random_state=1)

        #print('testing different models:')
        #models = []
        #models.append(('LR', LogisticRegression(solver='liblinear', multi_class='ovr')))
        #models.append(('LDA', LinearDiscriminantAnalysis()))
        #models.append(('KNN', KNeighborsClassifier()))
        #models.append(('CART', DecisionTreeClassifier(random_state=2)))  # secondBest
        #models.append(('NB', GaussianNB()))              #slow and bad
        #models.append(('SVM', SVC(gamma='auto'))) #slow
        #models.append(('RFC', RandomForestClassifier(random_state=1)))  # best

        #results = []
        #names = []
        #for name, model in models:
        #    kfold = StratifiedKFold(n_splits=10, random_state=1, shuffle=True)
        #    cv_results = cross_val_score(model, X_train, y_train, cv=kfold, scoring='accuracy')
        #    results.append(cv_results)
        #    names.append(name)
        #    print('%s: %f (%f)' % (name, cv_results.mean(), cv_results.std()))

        # DecisionTreeClassifier
        #   model = DecisionTreeClassifier()
        #   model.fit(X_train, y_train)
        #   predictions = model.predict(X_test)

        model = RandomForestClassifier(random_state=1)
        model.fit(X_train, y_train)
        predictions = model.predict(X_test)

        # Evaluate predictions
        # print(accuracy_score(y_test, predictions))
        print(confusion_matrix(y_test, predictions))
        #print(confusion_matrix(y_test, predictions))
        print(classification_report(y_test, predictions))
        print("Accuracy = " + str(round(accuracy_score(y_test, predictions) * 100, 2)) + " %")

        Plotter.plot_feature_importance_for_class(model, X_train)
        Plotter.plot_confusion_matrix(y_test, predictions)
        Plotter.plot_prediction_desicion(model, X_test, predictions, 0)
