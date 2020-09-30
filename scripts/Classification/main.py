from Classifier import Classifier
from DataManager import DataManager
import glob


def main():
    datamanager = DataManager()
    labels, data = datamanager.get_labels_and_data()
    classfier = Classifier(labels, data)
    classfier.predict()


if __name__ == "__main__":
    main()
