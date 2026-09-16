#!/usr/bin/env python3
"""
TSV Data Analyzer with Grouping, Percentage Calculation, and Visualization
Computes medians per group and generates barplots with standard deviation error bars
"""

__version__ = "1.0.0"

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
from typing import Optional, Dict, List
import sys


class TSVAnalyzer:
    """Analyze TSV data with optional grouping and percentage calculations"""
    
    def __init__(self, filepath: str):
        """
        Initialize the analyzer with a TSV file
        
        Args:
            filepath: Path to the TSV file
        """
        try:
            self.df = pd.read_csv(filepath, sep='\t')
            print(f"✓ Loaded TSV file with {len(self.df)} rows and {len(self.df.columns)} columns")
            print(f"  Columns: {', '.join(self.df.columns.tolist())}\n")
        except FileNotFoundError:
            print(f"Error: File '{filepath}' not found")
            sys.exit(1)
        except Exception as e:
            print(f"Error reading file: {e}")
            sys.exit(1)
    
    def calculate_percentages(self, value_columns: List[str], total_column: str) -> pd.DataFrame:
        """
        Calculate percentages based on a total column for multiple value columns
        Percentages are calculated BEFORE aggregation.
        
        Args:
            value_columns: List of column names to calculate percentages for
            total_column: Column to use as the denominator (total counts)
        
        Returns:
            DataFrame with value columns replaced by their percentage equivalents
        """
        # Validate all columns exist
        missing_cols = [col for col in value_columns + [total_column] if col not in self.df.columns]
        if missing_cols:
            print(f"Error: Columns not found: {missing_cols}")
            print(f"Available columns: {self.df.columns.tolist()}")
            sys.exit(1)
 
        df_copy = self.df.copy()
        
        print(f"✓ Calculated percentages relative to '{total_column}' (before aggregation):")
        for value_column in value_columns:
            df_copy[value_column] = (self.df[value_column] / self.df[total_column] * 100).round(2)
            print(f"  • '{value_column}': {df_copy[value_column].min():.2f}% - {df_copy[value_column].max():.2f}%")
        
        return df_copy
    
    def calculate_medians(self, group_column: str, value_columns: List[str]) -> pd.DataFrame:
        """
        Calculate medians per group for multiple value columns.
        Returns data in long format suitable for grouped barplot.
        
        Args:
            group_column: Column to group by
            value_columns: List of columns to calculate median for
        
        Returns:
            DataFrame in long format with medians and standard deviations per group and value
        """
        # Validate columns exist
        missing_cols = [col for col in value_columns if col not in self.df.columns]
        if missing_cols:
            print(f"Error: Columns not found: {missing_cols}")
            print(f"Available columns: {self.df.columns.tolist()}")
            sys.exit(1)
        
        if group_column not in self.df.columns:
            print(f"Error: Group column '{group_column}' not found")
            print(f"Available columns: {self.df.columns.tolist()}")
            sys.exit(1)
        
        all_stats = []
        
        for value_column in value_columns:
            stats = self.df.groupby(group_column)[value_column].agg([
                ('median', 'median'),
                ('std', 'std'),
                ('count', 'count'),
                ('mean', 'mean')
            ]).reset_index()
            
            stats['value'] = value_column
            all_stats.append(stats)
        
        # Combine all statistics
        combined_stats = pd.concat(all_stats, ignore_index=True)
        
        # Replace NaN std (from single values) with 0
        combined_stats['std'] = combined_stats['std'].fillna(0)
        
        # Reorder columns
        combined_stats = combined_stats[[group_column, 'value', 'median', 'std', 'count', 'mean']]
        
        print(f"✓ Calculated statistics for {len(value_columns)} value column(s) grouped by '{group_column}':")
        for vc in value_columns:
            subset = combined_stats[combined_stats['value'] == vc]
            print(f"\n  {vc}:")
            print(subset.to_string(index=False))
        print()
        
        return combined_stats
    
    def plot_barplot(self, 
                     stats_df: pd.DataFrame, 
                     group_column: str,
                     color_map: Optional[Dict[str, str]] = None,
                     output_file: Optional[str] = None,
                     title: Optional[str] = None):
        """
        Create a grouped barplot with error bars showing standard deviation.
        X-axis contains the value column names, hue variable is the group categories.
        
        Args:
            stats_df: DataFrame with calculated statistics (from calculate_medians)
            group_column: Name of the grouping column
            color_map: Dictionary mapping group values to colors (e.g., {'group_a': '#FF5733'})
            output_file: Path to save the figure (optional)
            title: Custom title for the plot (optional)
        """
        # Create figure with nice styling
        fig, ax = plt.subplots(figsize=(14, 6))
        
        value_cols = stats_df['value'].unique()
        groups = stats_df[group_column].unique()
        
        # Default color palette if not provided
        if color_map is None:
            colors = plt.cm.Set2(np.linspace(0, 1, len(groups)))
            color_map = {group: colors[i] for i, group in enumerate(groups)}
        
        # Ensure all groups have colors
        colors_list = [color_map.get(group, '#1f77b4') for group in groups]
        
        # Set up bar positions
        bar_width = 0.8 / len(groups)
        x_positions = np.arange(len(value_cols))
        
        # Plot bars for each group
        for idx, group in enumerate(groups):
            group_data = stats_df[stats_df[group_column] == group].copy()
            group_data = group_data.set_index('value').reindex(value_cols)
            
            medians = group_data['median'].fillna(0).values  # Fill NaN with 0 for display
            stds = group_data['std'].fillna(0).values
            
            x_pos = x_positions + (idx - len(groups)/2 + 0.5) * bar_width
            
            ax.bar(x_pos, medians, bar_width, 
                   label=group,
                   yerr=stds, capsize=5,
                   color=colors_list[idx], alpha=0.8, 
                   edgecolor='black', linewidth=1.2,
                   error_kw={'elinewidth': 1.5, 'ecolor': 'black', 'alpha': 0.7})
            
            # Add value labels on bars
            #for i, (median, std) in enumerate(zip(medians, stds)):
            #    if pd.notna(median):  # Use pandas notna for better NaN detection
            #        height = median
            #        ax.text(x_pos[i], height, f'{median:.1f}%', 
            #               ha='center', va='bottom', fontsize=8, fontweight='bold')
        
        # Customize plot
        ax.set_xlabel('', fontsize=12, fontweight='bold')
        ax.set_ylabel('Median percentage (+/- SD)', fontsize=12, fontweight='bold')
        ax.set_title(title or f'Median Values by {group_column}', 
                     fontsize=14, fontweight='bold', pad=20)
        ax.set_xticks(x_positions)
        ax.set_xticklabels(value_cols, fontsize=11)
        ax.legend(title=group_column, fontsize=10, title_fontsize=11, loc='best')
        ax.grid(axis='y', alpha=0.3, linestyle='--')
        
        plt.tight_layout()
        
        # Save if requested
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"✓ Plot saved to '{output_file}'")
        
        plt.show()
        return fig, ax


def parse_colors(color_string: str) -> Dict[str, str]:
    """
    Parse color string in format: 'group1:#FF5733,group2:#33FF57'
    
    Args:
        color_string: Comma-separated list of group:color pairs
    
    Returns:
        Dictionary mapping groups to colors
    """
    color_map = {}
    try:
        pairs = color_string.split(',')
        for pair in pairs:
            group, color = pair.strip().split(':')
            color_map[group.strip()] = color.strip()
        return color_map
    except ValueError:
        print("Error: Invalid color format. Use: 'group1:#FF5733,group2:#33FF57'")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description='Analyze TSV data with grouping, percentages, and visualization',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('tsv_file', help='Path to TSV file')
    parser.add_argument('--group', required=True, help='Column name for grouping variable')
    parser.add_argument('--value', action='append', required=True, dest='values',
                       help='Column name(s) for values to analyze (can use multiple times: --value col1 --value col2)')
    parser.add_argument('--percent', action='store_true', 
                       help='Calculate percentages (requires --total)')
    parser.add_argument('--total', help='Column name for total counts (used with --percent)')
    parser.add_argument('--colors', help='Custom colors as "group1:#HEX,group2:#HEX"')
    parser.add_argument('--output', help='Path to save the plot image')
    parser.add_argument('--title', help='Custom title for the plot')
    
    args = parser.parse_args()
    
    # Validation
    if args.percent and not args.total:
        print("Error: --total is required when using --percent")
        sys.exit(1)
    
    # Initialize analyzer
    analyzer = TSVAnalyzer(args.tsv_file)
    
    # Calculate percentages if requested (BEFORE aggregation)
    if args.percent:
        analyzer.df = analyzer.calculate_percentages(args.values, args.total)
        print(analyzer.df)
    # Calculate medians per group for all value columns at once
    print(f"{'='*60}")
    print(f"Calculating medians for {len(args.values)} value column(s)")
    print(f"{'='*60}\n")
    
    stats = analyzer.calculate_medians(args.group, args.values)
    
    # Parse custom colors if provided
    color_map = None
    if args.colors:
        color_map = parse_colors(args.colors)
        print(f"✓ Using custom colors: {color_map}\n")
    
    # Create single grouped plot
    analyzer.plot_barplot(stats, args.group, 
                         color_map=color_map, 
                         output_file=args.output,
                         title=args.title)


if __name__ == '__main__':
    main()
