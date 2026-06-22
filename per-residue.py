import pandas as pd
import matplotlib
# Force the 'Agg' backend to prevent errors in Conda/Colab environments
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import os

def generate_decomposition_plot(csv_file, output_png):
    if not os.path.exists(csv_file):
        print(f"❌ Error: {csv_file} not found.")
        return

    print(f"Reading {csv_file}...")
    
    # 1. Load the data 
    # Skipping the first few lines of metadata to reach the header
    try:
        df = pd.read_csv(csv_file, skiprows=3)
        # If 'Frame #' is not the first column, we try skipping 4 lines instead
        if 'Frame #' not in df.columns:
            df = pd.read_csv(csv_file, skiprows=4)
    except Exception as e:
        print(f"❌ Error reading CSV: {e}")
        return

    # 2. Clean the data
    # Remove any non-numeric rows (like 'Receptor:' headers inside the file)
    df = df[pd.to_numeric(df['TOTAL'], errors='coerce').notnull()]
    
    # Convert energy columns to numeric
    energy_cols = ['van der Waals', 'Electrostatic', 'Polar Solvation', 'Non-Polar Solv.', 'TOTAL']
    for col in energy_cols:
        df[col] = pd.to_numeric(df[col])

    # 3. Calculate Statistics (Mean and Std Dev) across all frames
    summary = df.groupby('Residue')['TOTAL'].agg(['mean', 'std']).reset_index()
    
    # Sort by mean energy (from most favorable to least favorable)
    summary = summary.sort_values(by='mean')

    # 4. Plotting
    plt.figure(figsize=(20, 8))
    
    # Color coding: Green for favorable binding, Red for unfavorable
    colors = ['#27ae60' if x < 0 else '#c0392b' for x in summary['mean']]
    
    # Create bar chart with error bars representing Standard Deviation
    plt.bar(summary['Residue'], summary['mean'], 
            yerr=summary['std'], 
            color=colors, 
            edgecolor='black', 
            alpha=0.8,
            error_kw={'capsize': 3, 'elinewidth': 0.6, 'ecolor': '#2c3e50'})

    # Chart Styling
    plt.axhline(0, color='black', linewidth=1.2)
    plt.ylabel('Binding Energy Contribution ($\Delta G$ kcal/mol)', fontweight='bold', fontsize=12)
    plt.xlabel('Residues (Chain:Type:Number)', fontweight='bold', fontsize=12)
    plt.title('Per-Residue Energy Decomposition: Average over Trajectory', fontsize=16, pad=25)
    
    # Formatting X-axis for 98 residues
    plt.xticks(rotation=90, fontsize=9)
    plt.grid(axis='y', linestyle='--', alpha=0.4)
    plt.tight_layout()

    # 5. Save Results
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    summary.to_csv('residue_decomposition_summary.csv', index=False)
    
    print("-" * 40)
    print(f"✅ Plot saved: {output_png}")
    print(f"📊 Summary data saved: residue_decomposition_summary.csv")
    print(f"Top 5 Binders Identified:")
    print(summary.head(5)[['Residue', 'mean']])

if __name__ == "__main__":
    # Ensure this matches your uploaded filename
    target_csv = "FINAL_DECOMP_MMPBSA.csv" 
    output_image = "mmpbsa_decomposition_chart.png"
    
    generate_decomposition_plot(target_csv, output_image)