import torch
import time

num_matrices = 1000000

# Set the default device to GPU
if torch.cuda.is_available():
    torch.cuda.set_device(0)  # You can specify the GPU index if you have multiple GPUs

# Function to generate random symmetric 12x12 matrices on the GPU
def generate_symmetric_matrices(batch_size):
    A = torch.randn(batch_size, 12, 12, device='cuda', dtype=torch.float64) * 2000 - 1000  # Generate random batch of 12x12 matrices on GPU
    A = (A + A.permute(0, 2, 1)) / 2  # Make them symmetric
    return A

# Number of matrices to generate in each batch
batch_size = 1000000  # Reduce batch size to fit GPU memory if necessary
# Number of batches
num_batches = num_matrices // batch_size

# Measure the time for EVD only
eigenvectors_time = 0

for _ in range(num_batches):
    matrices_batch = generate_symmetric_matrices(batch_size)
    
    # Measure the time for EVD
    start_time = time.time()
    eigenvalues_batch, eigenvectors_batch = torch.linalg.eigh(matrices_batch)
    end_time = time.time()
    
    eigenvectors_time += end_time - start_time

# Calculate the total time taken for EVD
print(f"Time taken for EVD of {num_matrices} matrices in {num_batches} batches on GPU: {eigenvectors_time} seconds")
print(f"Average time: {eigenvectors_time / num_batches} seconds")
