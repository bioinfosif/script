import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

def read_triangular_matrix(filename):
    """
    Read a triangular matrix from a TSV file and reconstruct the full symmetric matrix.
    
    Args:
        filename (str): Path to the input file
    
    Returns:
        pd.DataFrame: Symmetric matrix with proper labels
    """
    # Check if file exists
    if not Path(filename).exists():
        raise FileNotFoundError(f"File {filename} not found")
    
    # Read file line by line
    with open(filename, "r") as f:
        lines = [line.strip().split('\t') for line in f.readlines()]
    
    # Filter out empty lines
    lines = [line for line in lines if line and line[0]]
    
    # Extract labels and values
    labels = [row[0] for row in lines]
    values = [row[1:] for row in lines]
    
    n = len(labels)
    matrix = np.full((n, n), np.nan)
    
    # Fill lower triangular matrix
    for i in range(n):
        for j in range(min(i, len(values[i]))):
            try:
                val = float(values[i][j])
                if not np.isnan(val):  # Only assign non-NaN values
                    matrix[i][j] = val
            except (IndexError, ValueError):
                continue  # Skip invalid entries
    
    # Make matrix symmetric
    for i in range(n):
        for j in range(i + 1, n):
            if not np.isnan(matrix[j][i]):
                matrix[i][j] = matrix[j][i]
    
    # Fill diagonal (100 for similarity, 0 for distance)
    np.fill_diagonal(matrix, 100.0)
    
    # Create DataFrame
    df = pd.DataFrame(matrix, index=labels, columns=labels)
    return df.round(2)

def create_heatmap(df, output_file="matrix_heatmap_Staph.png", title="Symmetric Matrix Heatmap"):
    """
    Create and save a heatmap visualization of the matrix.
    
    Args:
        df (pd.DataFrame): Input matrix
        output_file (str): Output filename
        title (str): Plot title
    """
    # Set up the plot
    plt.figure(figsize=(40, 35))
    sns.set(font_scale=0.6)

    
    # Create heatmap with better styling
    mask = df.isnull()  # Mask NaN values
    sns.heatmap(
        df, 
        annot=True, 
        fmt=".1f",  # Show 1 decimal place for cleaner look
        cmap="RdYlBu_r",  # Better colormap for similarity matrices
        linewidths=0.5, 
        linecolor='white',
        mask=mask,
        square=True,  # Make cells square
        cbar_kws={"shrink": 0.6}
    )
    
    # Improve styling
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.title(title, fontsize=14, fontweight='bold', pad=20)
    
    # Add colorbar label
    cbar = plt.gcf().axes[-1]
    cbar.set_ylabel('Similarity Score', rotation=90, labelpad=15)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    
    return df

def analyze_matrix(df):
    """
    Provide basic statistics about the matrix.
    """
    print(f"Matrix dimensions: {df.shape}")
    print(f"Number of valid entries: {df.notna().sum().sum()}")
    print(f"Number of missing entries: {df.isna().sum().sum()}")
    
    # Exclude diagonal for statistics
    off_diagonal = df.values[np.triu_indices_from(df.values, k=1)]
    off_diagonal = off_diagonal[~np.isnan(off_diagonal)]
    
    if len(off_diagonal) > 0:
        print(f"\nOff-diagonal statistics:")
        print(f"Mean: {np.mean(off_diagonal):.2f}")
        print(f"Std: {np.std(off_diagonal):.2f}")
        print(f"Min: {np.min(off_diagonal):.2f}")
        print(f"Max: {np.max(off_diagonal):.2f}")

# Main execution
if __name__ == "__main__":
    filename = "Staph.csv"
    
    try:
        # Read and process the matrix
        print("Reading matrix...")
        df = read_triangular_matrix(filename)
        
        # Analyze the data
        print("\nMatrix Analysis:")
        analyze_matrix(df)
        
        # Create visualization
        print("\nCreating heatmap...")
        create_heatmap(
            df, 
            output_file="putida_similarity_heatmap.png",
            title="Putida Similarity Matrix"
        )
        
        print(f"\nHeatmap saved successfully!")
        print(f"Matrix shape: {df.shape}")
        
        # Optionally save the processed matrix
        df.to_csv("putida_processed_matrix.csv")
        print("Processed matrix saved as CSV")
        
    except Exception as e:
        print(f"Error: {e}")
        print("Please check your input file format and path.")
