import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("C://Users//gs010//Desktop//Bath//Year 4//PH40065 Final Year Project//Code//Path Integral//csv5.csv")

def plot_histogram(data, column, bins, title="", xlabel="", ylabel="Frequency"):
    plt.figure(figsize=(10, 6))
    plt.hist(data[column], bins=bins, color='blue', alpha=0.7, edgecolor='black')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(axis='y', alpha=0.75)
    plt.show()

plot_histogram(data, column="r", bins=200, title="Histogram to test randomness", xlabel="Random number")