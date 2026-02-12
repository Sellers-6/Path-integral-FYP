import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

L = 100

t = np.linspace(0, L, L)

def correlator(t):
    return 0.5 * np.exp(-t) + 0.5 * np.exp(t - (L - 1)) 

G = np.zeros(L)

for i in range (0, L):
    G[i] = correlator(i)

plt.plot(t, G)
plt.show()

# L = 500
# lambd = 1
# a = 0.5

# x = np.linspace(-1, 1, 2 * L + 1)

# def DPW(x):
#     return (lambd / 12) * (x**2 - a**2)**2

# V = DPW(x)

# fig, ax = plt.subplots()

# ax.plot(x, V, color="royalblue")

# # Hide top/right spines, move left/bottom to origin
# for side in ["top", "right"]:
#     ax.spines[side].set_visible(False)

# ax.spines["bottom"].set_position("zero")
# ax.spines["left"].set_position("zero")
# ax.spines["bottom"].set_color("black")
# ax.spines["left"].set_color("black")
# ax.spines["bottom"].set_linewidth(1.2)
# ax.spines["left"].set_linewidth(1.2)

# # Remove only the y-tick at 0 to avoid label overlap
# yticks = ax.get_yticks()
# ax.set_yticks([y for y in yticks if y != 0])

# # Label ±a
# ax.text(-a, -0.02, r"$-a$", ha="center", va="top", fontsize=12, color="black")
# ax.text(a, -0.02, r"$a$", ha="center", va="top", fontsize=12, color="black")

# # Manually reposition the x-axis label so it's not under the y-axis
# ax.set_xlabel("x", labelpad=10)
# ax.set_ylabel("V(x)", labelpad=10)
# ax.xaxis.set_label_coords(0.98, 0.23)   # move 'x' label slightly to the right

# # Title
# ax.set_title("Double-Well Potential", pad=10)

# plt.show()



# filename = "C:\\Users\gs010\Desktop\Bath\Year 4\PH40065 Final Year Project\Code\Path Integral\csv7.csv"  

# # Read the CSV
# data = pd.read_csv(filename)

# # Check that the data loaded correctly
# print(data.head())

# # Plot x(t)
# plt.plot(data["Time"], data["Position"], color="royalblue", linewidth=1.5)
# plt.xlabel("t")
# plt.ylabel("x(t)")
# plt.title("One path x(t) discretized by 50 time slices")
# plt.grid(True)
# plt.show()
