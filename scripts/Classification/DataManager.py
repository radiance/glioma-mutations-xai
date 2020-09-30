import pandas as pd
import numpy as np
from sklearn import preprocessing
from Preprocessor import Preprocessor

class DataManager:
    label_col = 'Age'

    data_frame = pd.DataFrame()
    data = pd.DataFrame()
    labels = pd.DataFrame()

    def __init__(self):
        self.__load_file()
        self.__process_data()
        self.__split_labels_from_data()

    def __load_file(self):
        dim_list = ['Age', 'TP53', 'IDH1', 'ATRX', 'NF1', 'PIK3R2', 'TERT', 'KMT2A', 'ERBB2', 'ATR', 'BCORL1', 'FLG', 'PIK3CG', 'KDM6A',
                    'IL7R', 'RAMP2', 'AXL', 'BRCA1', 'BARD1', 'PBRM1', 'EP300', 'KMT2D', 'RPTOR', 'U2AF1', 'EZH2', 'RELN',
                    'PHLPP1', 'SMO', 'CREBBP', 'KIT', 'H3F3B', 'USH2A', 'GNAS', 'SUZ12', 'MAP2K2', 'EPHA3', 'IRF4',
                    'SMARCA4', 'PTCH1', 'RYR2', 'EPHA5', 'PDGFRB', 'RAD21', 'CARD11', 'PPM1D', 'IRS2', 'SF3A1', 'SOX1',
                    'BLM', 'CDKN2A', 'FLT3', 'TCHH', 'SMC3', 'NTRK1', 'SUSD2', 'FUBP1', 'MED12', 'ASXL1', 'MUC4', 'PRDM1',
                    'FGFR3', 'PALB2', 'TEX13D', 'BCL6', 'APC', 'SVIL', 'ASXL2', 'ERBB3', 'MUTYH', 'TNFRSF14', 'OBSCN',
                    'RET', 'TNFAIP3', 'HMCN1', 'RICTOR', 'PCLO', 'FOXL2', 'ISM2', 'ABL1', 'CTNNB1', 'FGFR1', 'MPL', 'NF2',
                    'TSHR', 'ACVR1', 'MKI67', 'GRIN2A', 'JAK3']
        self.data_frame = pd.read_csv("..\data\mutations_merged_filtered_and_processed.csv", usecols=dim_list, sep=';')

    def __set_age_groups_as_label(self):
        for i, row in self.data_frame.iterrows():
            if self.data_frame.loc[i][self.label_col] <= 18:
                self.data_frame.at[i, self.label_col] = 1 #"0 >= age <= 23"
            elif self.data_frame.loc[i][self.label_col] <= 70:
                self.data_frame.at[i, self.label_col] = 2 #"23 > age <= 50"
            elif self.data_frame.loc[i][self.label_col] > 70:
                self.data_frame.at[i, self.label_col] = 3 #"50 > age"

    def __process_data(self):
        # # replace nan values with mean of column
        self.data_frame = Preprocessor.replaceNanValuesWithMedian(data_frame=self.data_frame)
        #self.data_frame = Preprocessor.countMutations(data_frame= self.data_frame)
        # # remove columns where majority of values is Nan or Zero
        # self.data_frame = Preprocessor.deleteNanColumns(data_frame=self.data_frame, threshold=99)

        #self.__set_text_as_true_and_nan_as_false()

        self.__set_age_groups_as_label()
        self.__categorical_data_to_numerical()

    def __categorical_data_to_numerical(self):
        dataframe_copy = self.data_frame.select_dtypes(include=['object']).copy()

        # join cols with int and float values
        int_dataframe = self.data_frame.select_dtypes(include=['int64']).copy()
        float_dataframe = self.data_frame.select_dtypes(include=['float64']).copy()
        dataframe_int_float = pd.concat([float_dataframe, int_dataframe], axis=1)

        le = preprocessing.LabelEncoder()
        dataframe_categorical = dataframe_copy.astype(str).apply(le.fit_transform)

        self.data_frame = pd.concat([dataframe_int_float, dataframe_categorical], axis=1)

    def get_whole_dataframe(self):
        return self.data_frame

    def __split_labels_from_data(self):
        dataframe = self.data_frame.copy()
        self.labels = dataframe.pop(self.label_col)
        self.data = dataframe

    def get_labels_and_data(self):
        return self.labels, self.data

    def __set_text_as_true_and_nan_as_false(self):
        additionalData = self.data_frame[["Age", 'Mutation_Count']].copy() #['Gender', "Age", 'Mutation_Count', 'DifferentMutatedGenesCount']
        self.data_frame = self.data_frame.drop(columns=["Age", 'DifferentMutatedGenesCount']) #["Age", 'Mutation_Count', 'DifferentMutatedGenesCount']
        self.data_frame.replace('0', np.nan, inplace=True)
        self.data_frame = pd.DataFrame(np.where(self.data_frame.isna(), self.data_frame, 0), columns=self.data_frame.columns)
        self.data_frame = self.data_frame.fillna(0)
        self.data_frame = self.data_frame.join(additionalData)
        self.data_frame.replace(np.nan, 0, inplace=True)
