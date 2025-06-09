# Does the visualizations of datasets be presented in the paper?
# Examine the original dataset
# Rewrite the description for each dataset
# Read the metrics, begin to put the execution metrics into paper

import pandas as pd
import argparse

import seaborn as sns
import matplotlib.pyplot as plt

from pathlib import Path



def main():
    
    # parser = argparse.ArgumentParser(description='Read CSV file')
    # parser.add_argument('args', nargs='+', help = 'CSV name')

    # args = parser.parse_args()

    folder_path = Path('../../data/UCR-time-series')
    file_list = [f.name[:-4] for f in folder_path.iterdir() if f.is_file()]
    for i in file_list:
        print(i)
        df = pd.read_csv(f'../../data/UCR-time-series/{i}.csv', index_col=0, header = None)
        df.columns = ['Value']
        df.index.name = 'No.'
        print(df.describe())
        print(df.var())

        plt.figure(figsize=(8, 4))
        sns.histplot(df['Value'], kde = True, bins=30, color='skyblue', edgecolor='black')
        plt.title('Histogram + KDE of Value')
        plt.xlabel('Value')
        plt.ylabel('Frequency')
        plt.grid(True)
        plt.savefig(f'visualization/{i}.png')
        plt.close()
        


if __name__ == '__main__':
    pd.set_option('display.float_format', '{:.6f}'.format)
    main()