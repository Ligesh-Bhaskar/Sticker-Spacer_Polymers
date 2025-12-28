import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import os
import glob
from scipy.signal import find_peaks
from scipy.integrate import trapezoid as trapz
import warnings
warnings.filterwarnings('ignore')
from matplotlib.patches import Patch

MONOMERS_PER_POLYMER = 12
NUMBER_OF_POLYMERS = 1000
POLYMER_RADIUS = 1.0   
SOLVENT_RADIUS = 4.0   
POLYMER_VOLUME = (4/3) * np.pi * POLYMER_RADIUS**3
SOLVENT_VOLUME = (4/3) * np.pi * SOLVENT_RADIUS**3

plot_params = {
    'label_fontsize': 40,
    'title_fontsize': 20,
    'legend_fontsize': 40,
    'tick_labelsize': 40,
    'border_thickness': 3,
    'major_tick_width': 3,
    'major_tick_length': 8,
    'minor_tick_width': 2,
    'minor_tick_length': 4,
    'major_tick_pad': 10,
    'dpi': 300
}

plt.rcParams.update({
    'figure.figsize': (12, 8),
    'font.size': plot_params['tick_labelsize'],
    'axes.labelsize': plot_params['label_fontsize'],
    'axes.titlesize': plot_params['title_fontsize'],
    'xtick.labelsize': plot_params['tick_labelsize'],
    'ytick.labelsize': plot_params['tick_labelsize'],
    'legend.fontsize': plot_params['legend_fontsize'],
    'axes.linewidth': plot_params['border_thickness'],
    'xtick.major.width': plot_params['major_tick_width'],
    'xtick.major.size': plot_params['major_tick_length'],
    'ytick.major.width': plot_params['major_tick_width'],
    'ytick.major.size': plot_params['major_tick_length'],
    'xtick.minor.width': plot_params['minor_tick_width'],
    'xtick.minor.size': plot_params['minor_tick_length'],
    'ytick.minor.width': plot_params['minor_tick_width'],
    'ytick.minor.size': plot_params['minor_tick_length'],
    'xtick.major.pad': plot_params['major_tick_pad'],
    'ytick.major.pad': plot_params['major_tick_pad'],
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'font.family': 'sans-serif',
    'font.weight': 'bold'    
})

class CombinedAnalyzer:
    
    def __init__(self, base_dir='PostProcessing', pe_values=[0, 0.5, 1.0, 2.0]):
		
        self.base_dir = base_dir
        self.cluster_dir = os.path.join(base_dir, 'Cluster_Data')
        self.partition_dir = os.path.join(base_dir, 'Partition_Data')
        self.rg_dir = os.path.join(base_dir, 'Rg_Data')
        self.output_dir = os.path.join(base_dir, 'Plots')
        self.pe_values = pe_values
        self.colors = plt.cm.viridis(np.linspace(0, 1, len(pe_values)))
        self.plot_params = plot_params
        os.makedirs(self.output_dir, exist_ok=True)
        
    def ParseFilename(self, filename):
		
        params = {}
        parts = filename.split('_')
        
        for part in parts:
            if part.startswith('T'):
                params['T'] = int(part[1:])
            elif part.startswith('ABP'):
                params['ABP'] = int(part[3:])
            elif part.startswith('Pe'):
                params['Pe'] = float(part[2:])
            elif part.startswith('E'):
                e_part = part[1:]
                if e_part.replace('.', '').replace('-', '').isdigit():
                    params['E'] = float(e_part)
        
        return params
    
    def ConfigureAxes(self, ax):
		
        ax.tick_params(which='major', width=self.plot_params['major_tick_width'], 
                      length=self.plot_params['major_tick_length'], direction='in')
        ax.tick_params(which='minor', width=self.plot_params['minor_tick_width'], 
                      length=self.plot_params['minor_tick_length'], direction='in')
        for spine in ax.spines.values():
            spine.set_linewidth(self.plot_params['border_thickness'])
    
    def LoadClusterStatistics(self, t_value, e_value, abp_value, pe_value):
		
        pattern = f"Cluster_Statistics_T{t_value}_ABP{abp_value}_Pe{pe_value}_E{e_value}_*.dat"
        files = glob.glob(os.path.join(self.cluster_dir, pattern))
        
        if not files:
            print(f"Warning: No cluster file found for T{t_value}, E{e_value}, ABP{abp_value}, Pe{pe_value}")
            return None
            
        filename = files[0] 
        try:
            data = np.loadtxt(filename, skiprows=1)
            if data.ndim == 1:
                data = data.reshape(1, -1)
            return {
                'frame': data[:, 0],
                'avg_cluster_size': data[:, 1],
                'std_cluster_size': data[:, 2],
                'avg_polymer_cluster_size': data[:, 3],
                'std_polymer_cluster_size': data[:, 4],
                'largest_cluster': data[:, 5],
                'num_clusters': data[:, 6],
                'num_polymer_clusters': data[:, 7],
                'num_mixed_clusters': data[:, 8],
                'num_solvent_clusters': data[:, 9],
                'fraction_largest': data[:, 10],
                'weight_avg_cluster_size': data[:, 11]
            }
        except Exception as e:
            print(f"Error loading {filename}: {e}")
            return None
    
    def ReadRgData(self, filepath):
		
        try:
            data = pd.read_csv(filepath, sep=r'\s+', comment='#',
                             names=['Dense_Phase_Rg', 'Dense_Asphericity', 'Dense_Acylindricity', 
                                   'Dense_RelShapeAniso', 'Dense_Lambda2_Lambda1', 'Dense_Lambda3_Lambda1',
                                   'Dilute_Phase_Rg', 'Dilute_Asphericity', 'Dilute_Acylindricity', 
                                   'Dilute_RelShapeAniso', 'Dilute_Lambda2_Lambda1', 'Dilute_Lambda3_Lambda1',
                                   'Largest_Cluster_Rg', 'Largest_Asphericity', 'Largest_Acylindricity', 
                                   'Largest_RelShapeAniso', 'Largest_Lambda2_Lambda1', 'Largest_Lambda3_Lambda1'])
            return data
        except Exception as e:
            print(f"Error reading Rg data {filepath}: {e}")
            return None
    
    def ReadPartitionData(self, filepath):
		
        try:
            if 'Summary' in filepath:
                data = pd.read_csv(filepath, sep=r'\s+', comment='#', skiprows=3,
                                 names=['Parameter', 'Average_Value', 'Std_Dev'])
                return data
            else:
                data = pd.read_csv(filepath, sep=r'\s+', comment='#',
                                 names=['Frame_Num', 'Polymer_K', 'Solvent_K', 'Mixture_K',
                                       'Dense_Polymers', 'Dilute_Polymers', 
                                       'Dense_Solvents', 'Dilute_Solvents', 'Mixture_K_Weighted'])                                     
                return data
        except Exception as e:
            print(f"Error reading partition data {filepath}: {e}")
            return None
    
    def CollectPartitionData(self):
		
        partition_data = []
        scale_fact = 0.25        
        partition_files = glob.glob(os.path.join(self.partition_dir, '*.dat'))
        
        for filepath in partition_files:
            filename = os.path.basename(filepath)
            params = self.ParseFilename(filename)
            
            if params and 'Summary' not in filename:
                part_data = self.ReadPartitionData(filepath)
                if part_data is not None:
                    polymer_k = part_data['Polymer_K'].mean()
                    polymer_k_std = part_data['Polymer_K'].std()*scale_fact
                    solvent_k = part_data['Solvent_K'].mean()
                    solvent_k_std = part_data['Solvent_K'].std()*scale_fact
                    mixture_k = part_data['Mixture_K'].mean()
                    mixture_k_std = part_data['Mixture_K'].std()*scale_fact
                    
                    result = {**params, 'Polymer_K': polymer_k, 'Polymer_K_std': polymer_k_std,
                             'Solvent_K': solvent_k, 'Solvent_K_std': solvent_k_std,
                             'Mixture_K': mixture_k, 'Mixture_K_std': mixture_k_std,
                             'filepath': filepath}
                    partition_data.append(result)
        
        return pd.DataFrame(partition_data)
    
    def PlotNormalizedAvgClusterSize(self, t_values, e_values, abp_values, norm_value=2200, top_n=8):
		
        pe_colors = {0: 'black', 0.5: 'olive', 1.0: 'blue', 2.0: 'red'}
        
        for t_val in t_values:
            for e_val in e_values:
                for abp_val in abp_values:
                    fig, ax = plt.subplots(figsize=(15, 10))
                    
                    pe_means = []
                    pe_stds = []
                    pe_labels = []
                    bar_colors = []
                    
                    for i, pe_val in enumerate(self.pe_values):
                        cluster_data = self.LoadClusterStatistics(t_val, e_val, abp_val, pe_val)
                        if cluster_data is None:
                            continue
                        
                        avg_topN_per_frame = []
                        
                        for frame_idx in range(len(cluster_data['frame'])):
                            largest_cluster = cluster_data['largest_cluster'][frame_idx]
                            num_clusters = cluster_data['num_clusters'][frame_idx]
                            weight_avg = cluster_data['weight_avg_cluster_size'][frame_idx]
                    
                            if num_clusters >= top_n:
                                estimated_topN_avg = min(largest_cluster, weight_avg * 1.5)
                            else:
                                estimated_topN_avg = weight_avg
                            
                            avg_topN_per_frame.append(estimated_topN_avg)
                        
                        if not avg_topN_per_frame:
                            continue
                        
                        mean_topN = np.mean(avg_topN_per_frame)
                        std_topN = np.std(avg_topN_per_frame)
                        
                        normalized_mean = mean_topN / norm_value
                        normalized_std = std_topN / norm_value
                        
                        pe_means.append(normalized_mean)
                        pe_stds.append(normalized_std)
                        pe_labels.append(str(pe_val))
                        bar_colors.append(pe_colors.get(pe_val, 'gray'))
                    
                    if not pe_means:  
                        plt.close()
                        continue
                    
                    bars = ax.bar(pe_labels, pe_means, yerr=pe_stds, color=bar_colors, 
                                  alpha=0.9, capsize=10, error_kw={'elinewidth': 3.0, 'capthick': 3.0})
                    
                    ax.set_xlabel(r'$Pe$')
                    ax.set_ylabel(r'$Average \, Size$') 
                    self.ConfigureAxes(ax)
                    
                    plt.tight_layout()
                    filename = f'{self.output_dir}/Normalized_avg_clustersize_T{t_val}_E{e_val}_ABP{abp_val}_Pe_comparison.png'
                    plt.savefig(filename, dpi=self.plot_params['dpi'], bbox_inches='tight')
                    plt.close()
    
    def PlotPartitionCoefficientsParameters(self, partition_df, pe_values, epsilon_values, abp_values):
		
        for abp_val in abp_values:
            abp_data = partition_df[partition_df['ABP'] == abp_val]
            
            for epsilon_val in epsilon_values:
                epsilon_data = abp_data[abp_data['E'] == epsilon_val]
                
                if epsilon_data.empty:
                    continue
                
                plot_data = []
                
                for pe_val in pe_values:
                    pe_data = epsilon_data[epsilon_data['Pe'] == pe_val]
                    
                    if not pe_data.empty:
                        polymer_means = pe_data['Polymer_K'].values
                        polymer_stds = pe_data['Polymer_K_std'].values
                        
                        solvent_means = pe_data['Solvent_K'].values
                        solvent_stds = pe_data['Solvent_K_std'].values
                        
                        mixture_means = pe_data['Mixture_K'].values
                        mixture_stds = pe_data['Mixture_K_std'].values
                        
                        if len(polymer_means) == 1:
                            polymer_final_mean = polymer_means[0]
                            polymer_final_std = polymer_stds[0]
                            
                            solvent_final_mean = solvent_means[0]
                            solvent_final_std = solvent_stds[0]
                            
                            mixture_final_mean = mixture_means[0]
                            mixture_final_std = mixture_stds[0]
                        else:
                            polymer_final_mean = np.mean(polymer_means)
                            polymer_final_std = np.sqrt(np.sum(polymer_stds**2)) / len(polymer_stds)
                            
                            solvent_final_mean = np.mean(solvent_means)
                            solvent_final_std = np.sqrt(np.sum(solvent_stds**2)) / len(solvent_stds)
                            
                            mixture_final_mean = np.mean(mixture_means)
                            mixture_final_std = np.sqrt(np.sum(mixture_stds**2)) / len(mixture_stds)
                                                    
                        plot_data.append({
                            'Pe': pe_val,
                            'Polymer_K_mean': polymer_final_mean,
                            'Polymer_K_std': polymer_final_std,
                            'Solvent_K_mean': solvent_final_mean,
                            'Solvent_K_std': solvent_final_std,
                            'Mixture_K_mean': mixture_final_mean,
                            'Mixture_K_std': mixture_final_std
                        })
                
                if not plot_data:
                    continue
                    
                pe_groups = pd.DataFrame(plot_data)
                pe_groups = pe_groups.fillna(0)
                
                fig, ax = plt.subplots(figsize=(15, 10))
                
                pe_vals = pe_groups['Pe'].values
                x_pos = np.arange(len(pe_vals))
                width = 0.4
                
                max_value = max(pe_groups['Polymer_K_mean'].max() + pe_groups['Polymer_K_std'].max(),
                                pe_groups['Solvent_K_mean'].max() + pe_groups['Solvent_K_std'].max())
                ax.set_ylim(0, max_value * 1.2)
                
                bars2 = ax.bar(x_pos, pe_groups['Polymer_K_mean'], width, yerr=pe_groups['Polymer_K_std'], capsize=10, 
                               color='purple', alpha=0.9, label='Polymer', error_kw={'elinewidth': 3.0, 'capthick': 3.0})
                
                bars3 = ax.bar(x_pos + width, pe_groups['Solvent_K_mean'], width, yerr=pe_groups['Solvent_K_std'], capsize=10,
                               color='darkorange', alpha=0.9, label='Bath', error_kw={'elinewidth': 3.0, 'capthick': 3.0})               
                
                ax.set_xlabel(r'$Pe$')
                ax.set_ylabel(r'$Average \, Partition \, Coefficient$', fontsize=34)
                ax.set_xticks(x_pos)
                ax.set_xticklabels([f'{pe}' for pe in pe_vals])
                ax.legend(frameon=False)
                ax.grid(False)
                
                plt.tight_layout()
                plt.savefig(os.path.join(self.output_dir, f'Partition_coefficients_ABP{abp_val}_E{epsilon_val}.png'), 
                            dpi=self.plot_params['dpi'], bbox_inches='tight')
                plt.close()
    
    def PlotRgHistogramsParameters(self, pe_values, epsilon_values, abp_values):
		
        for abp_val in abp_values:
            for epsilon_val in epsilon_values:
                
                fig, axes = plt.subplots(2, 2, figsize=(20, 15))
                axes = axes.flatten()
                
                for i, pe_val in enumerate(pe_values):
                    if i >= 4:
                        break
                    
                    rg_pattern = f"Individual_Rg_T*_ABP{abp_val}_Pe{pe_val}_E{epsilon_val}_*.dat"
                    rg_files = glob.glob(os.path.join(self.rg_dir, rg_pattern))
                    
                    if not rg_files:
                        continue
                    
                    dense_rg_values = []
                    dilute_rg_values = []
                    
                    for rg_file in rg_files:
                        rg_data = self.ReadRgData(rg_file)
                        if rg_data is not None:
                            dense_rg = rg_data['Dense_Phase_Rg'].dropna()
                            dense_rg = dense_rg[dense_rg > 0]
                            
                            dilute_rg = rg_data['Dilute_Phase_Rg'].dropna()
                            dilute_rg = dilute_rg[dilute_rg > 0]
                            
                            dense_rg_values.extend(dense_rg)
                            dilute_rg_values.extend(dilute_rg)
                    
                    if dense_rg_values or dilute_rg_values:
                        if dense_rg_values:
                            axes[i].hist(dense_rg_values, alpha=0.7, bins=20, 
                                       color='magenta', label='Dense', density=True)
                        if dilute_rg_values:
                            axes[i].hist(dilute_rg_values, alpha=0.7, bins=20, 
                                       color='green', label='Dilute', density=True)
                        
                        axes[i].set_xlabel(r'$R_g$')
                        axes[i].set_xlim(right=12)
                        axes[i].set_ylabel(r'$P(R_g)$')
                        axes[i].legend(frameon=False)
                        axes[i].grid(False)
                        axes[i].text(0.95, 0.5, f'Pe = {pe_val}\n$\\epsilon$ = {epsilon_val}', transform=axes[i].transAxes, fontsize=self.plot_params['tick_labelsize'],
                                    verticalalignment='top', horizontalalignment='right', bbox=dict(facecolor='none', edgecolor='none'))
                        
                for i in range(len(pe_values), 4):
                    axes[i].set_visible(False)
                
                plt.tight_layout()
                plt.savefig(os.path.join(self.output_dir, f'Rg_Histograms_ABP{abp_val}_E{epsilon_val}.png'), 
                            dpi=self.plot_params['dpi'], bbox_inches='tight')
                plt.close()
    
    def RunAnalysisPlots(self, t_values, abp_values, pe_values, eps_values, pe_values_rg, eps_values_rg):
        
        print("\n[1/3] Generating Cluster Size Plots...")
        self.PlotNormalizedAvgClusterSize(t_values, eps_values, abp_values)
        
        print("\n[2/3] Generating Partition Coefficient Plots...")
        partition_df = self.CollectPartitionData()
        if not partition_df.empty:
            self.PlotPartitionCoefficientsParameters(partition_df, pe_values, eps_values, abp_values)
        else:
            print("Warning: No partition data found")
        
        print("\n[3/3] Generating Rg Histograms...")
        self.PlotRgHistogramsParameters(pe_values_rg, eps_values_rg, abp_values)
        
        print(f"\nAll Plots Saved to: {self.output_dir}")
        
        return partition_df

def main():

    t_values = [300]
    abp_values = [1200]    
    pe_values = [0, 0.5, 1.0, 2.0]
    eps_values = [6]
    pe_values_rg = [0, 0.5]     
    eps_values_rg = [4, 6]
    
    analyzer = CombinedAnalyzer(base_dir='.', pe_values=pe_values)
    
    analyzer.RunAnalysisPlots(t_values, abp_values, pe_values, eps_values, pe_values_rg, eps_values_rg)
    
if __name__ == "__main__":
    main()
