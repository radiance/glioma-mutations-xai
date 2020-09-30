import pandas as pd
import math


class Preprocessor:

    @staticmethod
    def isNanOrZero(value):
        if isinstance(value, str):
            if "nan" in value.lower() or "0" in value:
                return True
        elif isinstance(value, float):
            if math.isnan(value) or value == 0.0:
                return True
        elif isinstance(value, int):
            if value == 0:
                return True
        return False

    @staticmethod
    def deleteRowIfColumnIsNan(data_frame: pd.DataFrame, column_name: str):
        for index, row in data_frame.iterrows():
            if not (row[column_name] > 0):
                data_frame.drop(index, inplace=True)
        return data_frame

    @staticmethod
    def replaceNanValuesWithMedian(data_frame: pd.DataFrame):
        return data_frame.fillna(data_frame.median())

    @staticmethod
    def findAgeColumnName(data_frame: pd.DataFrame) -> str:
        for column_name in data_frame:
            if "age" in column_name.lower():
                return column_name

    @staticmethod
    def countMutations(data_frame: pd.DataFrame):
        for column_name in data_frame:
            if column_name not in ["Age", 'Mutation_Count']:
                data_frame[column_name] = data_frame[column_name].astype(str).str.split(' ').str.len()
        return data_frame

    # delete nan columns depending on the threshold how many (in percentage) of the values can be missing in column
    @staticmethod
    def deleteNanColumns(data_frame: pd.DataFrame, threshold: float = 0.0):
        missing_values_counter = 0
        columns_to_delete = list()
        for column_name in data_frame:
            column_size = 0
            missing_values_counter = 0
            for value in data_frame[column_name]:
                column_size += 1
                if Preprocessor.isNanOrZero(value):
                    missing_values_counter += 1
            percentage_of_missing_values = (missing_values_counter / column_size) * 100

            if percentage_of_missing_values > threshold:
                if column_name not in ["Age", 'Mutation_Count']:
                    columns_to_delete.append(column_name)
        for column_name in columns_to_delete:
            data_frame = data_frame.drop(column_name, 1)

        return data_frame
