time_txt_content = open("out.txt").read()
# Correcting the parsing logic
lines = time_txt_content.strip().split("\n")
# Resetting lists
matrix_sizes = []
times_eigen, errors_eigen = [], []
times_ours, errors_ours = [], []
times_cusolver, errors_cusolver = [], []
# Your colors data
colors = {
    "Ours": "#282F44",
    "Sympy": "#495D63",
    "Eigen": "#E28413",
    "Manual": "#C42847",
    "cuBlas": "#76B900"
}

i = 0
while i < len(lines):
    if lines[i].startswith("Matrix size:"):
        matrix_size = int(lines[i].split(":")[1].strip())
        matrix_sizes.append(matrix_size)
    elif "Eigen's eigen value decomposition" in lines[i]:
        times_eigen.append(float(lines[i].split(":")[1].split()[0]))
        errors_eigen.append(float(lines[i+1].split(":")[1]))
    elif "our eigen value decomposition" in lines[i]:
        times_ours.append(float(lines[i].split(":")[1].split()[0]))
        errors_ours.append(float(lines[i+1].split(":")[1]))
    elif "cusolver eigen value decomposition" in lines[i]:
        times_cusolver.append(float(lines[i].split(":")[1].split()[0]))
        errors_cusolver.append(float(lines[i+2].split(":")[1]))
        i += 1 # Skip the "Batched syevj execution successful!" line
    i += 1
import matplotlib.pyplot as plt
# Plot for the speed with integer X-axis
plt.figure(figsize=(8, 6))
plt.plot(matrix_sizes, times_eigen, color=colors["Eigen"], label="Eigen")
plt.plot(matrix_sizes, times_ours, color=colors["Ours"], label="Ours")
plt.plot(matrix_sizes, times_cusolver, color=colors["cuBlas"], label="cuSOLVER")
plt.xlabel("Matrix Size")
plt.ylabel("Execution Time (ms)")
# plt.title("Speed Comparison")
plt.legend()
# plt.grid(True)
plt.xticks(matrix_sizes)

plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig("evd_speed.png")
plt.close()

# Plot for the accuracy in log scale with multiple lines for each method and integer X-axis
plt.figure(figsize=(8, 6))
plt.semilogy(matrix_sizes, errors_eigen, color=colors["Eigen"], label="Eigen")
plt.semilogy(matrix_sizes, errors_ours, color=colors["Ours"], label="Ours")
plt.semilogy(matrix_sizes, errors_cusolver, color=colors["cuBlas"], label="cuSOLVER")
plt.xlabel("Matrix Size")
plt.ylabel("Maximum Error (log scale)")
# plt.title("Accuracy Comparison")
plt.legend()
# plt.grid(True)
plt.xticks(matrix_sizes)

plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig("evd_accuracy.png")
