import os
import pandas as pd
import matplotlib
# Force non-interactive backend to avoid ValueError in Conda/Colab
matplotlib.use('Agg') 
import matplotlib.pyplot as plt

# Function to extract frames from CSV files
def extract_frames(csv_files):
    frames = set()
    for csv_file in csv_files:
        df = pd.read_csv(csv_file)
        # Using a safer split/strip logic for PDB_X_ns naming
        frames.update(df['pdbid'].apply(lambda x: int(x.split('_')[1].replace('ns', ''))))
    return sorted(frames)

# Function to create timeline
def create_timeline(frames, csv_files):
    timeline = {}
    for frame in frames:
        timeline[frame] = {}
        for csv_file in csv_files:
            df = pd.read_csv(csv_file)
            # Match the specific ID format in your CSVs
            if f'PDB_{frame}_PROTEIN' in df['pdbid'].values:
                bond_type = os.path.splitext(os.path.basename(csv_file))[0]
                if bond_type not in timeline[frame]:
                    timeline[frame][bond_type] = 1
                else:
                    timeline[frame][bond_type] += 1
        if not timeline[frame]:
            timeline[frame] = {'No Bond': 1}
    return timeline

# Main function
def main():
    # Ensure this folder exists in your current directory
    bond_csv_folder = "bond_csv_files"
    if not os.path.exists(bond_csv_folder):
        print(f"Error: Folder '{bond_csv_folder}' not found.")
        return

    csv_files = [os.path.join(bond_csv_folder, file) for file in os.listdir(bond_csv_folder) if file.endswith('.csv')]

    frames = extract_frames(csv_files)
    timeline = create_timeline(frames, csv_files)

    # Prepare data for stacked bar chart
    data = {bond_type: [] for bond_type in csv_files}
    for frame in timeline:
        for bond_type in data:
            data[bond_type].append(timeline[frame].get(os.path.splitext(os.path.basename(bond_type))[0], 0))

    # Plot stacked bar chart
    plt.figure(figsize=(15, 6))
    bottom = None
    colors = {
        'hydrophobic_interactions': 'grey',
        'hydrogen_bonds': 'blue',
        'salt_bridges': 'yellow',
        'halogen_bonds': 'cyan',
        'pi_cation_interactions': 'maroon',
        'pi_stacks': 'green'
    }

    for bond_type in data:
        # Clean label for legend
        clean_label = bond_type.replace('_', ' ').title().split("/")[-1].split(".")[0]
        color_key = os.path.splitext(os.path.basename(bond_type))[0]
        
        if bottom is None:
            plt.bar(frames, data[bond_type], label=clean_label, color=colors.get(color_key, 'black'))
            bottom = data[bond_type]
        else:
            plt.bar(frames, data[bond_type], bottom=bottom, label=clean_label, color=colors.get(color_key, 'black'))
            bottom = [bottom[i] + data[bond_type][i] for i in range(len(bottom))]

    plt.xlabel('Time (ns)')
    plt.ylabel('Presence of Bonds')
    plt.title('Bond Timeline')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    # Adjusted ticks based on your 1000ns scale
    plt.xticks(range(0, 10, 1))
    plt.tight_layout()

    # SAVE TO CURRENT DIRECTORY
    output_filename = "bond_timeline_plot.png"
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Success: Plot saved as {output_filename} in {os.getcwd()}")
    plt.close()

if __name__ == "__main__":
    main()