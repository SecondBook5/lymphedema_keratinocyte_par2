"""
Publication-quality visualization for scRNA-seq analysis.
Standards: JCI Insight / Nature figure specifications.

Figure sizes respect JCI Insight limits (max 6.85 x 8.65 inches).
Fonts embedded as TrueType (PDF fonttype 42) for Illustrator editing.
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import seaborn as sns
from cycler import cycler
from pathlib import Path
from scipy import stats

from support import analysis_config as ac

#  Figure dimensions (inches) 
FIGSIZE = {
    'single':       (3.4, 3.2),
    'single_tall':  (3.4, 4.0),
    'double':       (6.85, 3.2),
    'double_tall':  (6.85, 5.5),
    'triple':       (6.85, 2.8),
    'umap':         (3.2, 3.0),
    'umap_pair':    (6.5, 3.0),
    'umap_2x2':     (6.5, 6.0),
    'dotplot':      (7.5, 4.5),
    'violin':       (3.4, 3.8),
    'violin_wide':  (5.5, 3.8),
    'volcano':      (3.8, 4.2),
    'heatmap':      (6.85, 5.0),
}

#  Color system 
COLORS = {
    # Primary condition colors
    ac.CASE_CONDITION:     ac.CONDITION_COLORS[ac.CASE_CONDITION],
    ac.CONTROL_CONDITION:  ac.CONDITION_COLORS[ac.CONTROL_CONDITION],
    'case':                ac.CONDITION_COLORS[ac.CASE_CONDITION],
    'control':             ac.CONDITION_COLORS[ac.CONTROL_CONDITION],

    # Cell type palette  perceptually distinct, colorblind-aware
    # Mast cells in distinctive orange: the key cell type of this paper
    'cell_types': {
        'Keratinocytes':        '#4E9A8F',   # Teal
        'Fibroblasts':          '#7B5EA7',   # Purple
        'Endothelial':          '#E8749A',   # Rose
        'Pericytes':            '#BDA0CB',   # Lavender
        'Melanocytes':          '#8B6347',   # Brown
        'Schwann cells':        '#9EB4C4',   # Steel blue
        'T cells':              '#D94F3D',   # Red
        'CD4+ T cells':         '#E07B54',   # Salmon
        'CD8+ T cells':         '#B03A2E',   # Dark red
        'Regulatory T cells':   '#F0A87C',   # Peach
        'NK cells':             '#E8A020',   # Amber
        'B cells':              '#3A7EC9',   # Blue
        'Plasma cells':         '#7FB8E8',   # Light blue
        'Macrophages':          '#2A7A30',   # Forest green
        'Monocytes':            '#5CB85C',   # Medium green
        'Dendritic cells':      '#9ABF5E',   # Yellow-green
        'Mast cells':           '#FF6B35',   # Bright orange  KEY CELL TYPE
        'Neutrophils':          '#FFB347',   # Light orange
        'Unknown':              '#AAAAAA',   # Gray
    },

    # 20-color palette for Leiden clusters
    'categorical': [
        '#D94F3D', '#3A7EC9', '#4E9A8F', '#7B5EA7', '#E8A020',
        '#2A7A30', '#E8749A', '#FF6B35', '#7FB8E8', '#9ABF5E',
        '#8B6347', '#BDA0CB', '#E07B54', '#5CB85C', '#9EB4C4',
        '#B03A2E', '#FFB347', '#4DB6AC', '#A66BBE', '#81C784',
    ],

    # Colormaps
    'expression':  'magma',    # Gene expression on UMAPs
    'diverging':   'RdBu_r',   # Fold change / centered data
    'sequential':  'YlOrRd',   # Counts / density
}

#  Typography 
FONTSIZE = {
    'panel_label':  11,    # Bold A, B, C labels
    'title':         9,
    'axis_title':    8,
    'tick':          7,
    'legend':        7,
    'annotation':    6.5,
}

#  Base rcParams 
_RCPARAMS = {
    'figure.dpi':             150,
    'savefig.dpi':            300,
    'font.family':            'sans-serif',
    'font.sans-serif':        ['Arial', 'Helvetica Neue', 'Helvetica', 'DejaVu Sans'],
    'font.size':              FONTSIZE['tick'],
    'axes.labelsize':         FONTSIZE['axis_title'],
    'axes.titlesize':         FONTSIZE['title'],
    'xtick.labelsize':        FONTSIZE['tick'],
    'ytick.labelsize':        FONTSIZE['tick'],
    'legend.fontsize':        FONTSIZE['legend'],
    'legend.title_fontsize':  FONTSIZE['legend'],
    'legend.frameon':         False,
    'legend.borderpad':       0.3,
    'axes.linewidth':         0.75,
    'axes.spines.top':        False,
    'axes.spines.right':      False,
    'axes.edgecolor':         '#2C2C2C',
    'axes.labelcolor':        '#2C2C2C',
    'text.color':             '#2C2C2C',
    'xtick.major.width':      0.75,
    'ytick.major.width':      0.75,
    'xtick.major.size':       3.0,
    'ytick.major.size':       3.0,
    'xtick.color':            '#2C2C2C',
    'ytick.color':            '#2C2C2C',
    'lines.linewidth':        1.25,
    'patch.linewidth':        0.5,
    'axes.grid':              False,
    # TrueType fonts embedded in PDF/SVG  editable in Illustrator/Inkscape
    'pdf.fonttype':           42,
    'ps.fonttype':            42,
    'svg.fonttype':           'none',
}


def set_publication_style():
    """Apply publication-quality styling globally."""
    mpl.rcParams.update(_RCPARAMS)
    plt.rcParams['axes.prop_cycle'] = cycler('color', COLORS['categorical'])
    sns.set_style('ticks', {
        'axes.linewidth':    0.75,
        'axes.edgecolor':    '#2C2C2C',
        'xtick.major.width': 0.75,
        'ytick.major.width': 0.75,
    })
    sns.set_context('paper', font_scale=1.0)


#  Color helpers 

def get_condition_colors():
    """Return condition  color mapping."""
    return dict(ac.CONDITION_COLORS)


def get_condition_order(values=None):
    """Return configured condition order, optionally filtered to observed values."""
    return ac.get_condition_order(values)


def get_case_condition():
    """Return the configured case condition label."""
    return ac.CASE_CONDITION


def get_control_condition():
    """Return the configured control condition label."""
    return ac.CONTROL_CONDITION


def get_contrast_label(joiner='vs'):
    """Return the configured human-readable contrast label."""
    return ac.contrast_label(joiner=joiner)


def get_delta_label():
    """Return the configured delta label."""
    return ac.delta_label()


def get_contrast_stem():
    """Return the configured filesystem-safe contrast stem."""
    return ac.contrast_stem()


def get_cell_type_colors(cell_types=None):
    """Return cell type  color mapping, assigning categorical fallbacks."""
    palette = dict(COLORS['cell_types'])
    if cell_types is None:
        return palette
    out = {}
    for i, ct in enumerate(cell_types):
        out[ct] = palette.get(ct, COLORS['categorical'][i % len(COLORS['categorical'])])
    return out


#  Panel labels 

def add_panel_label(ax, label, x=-0.15, y=1.06):
    """Add bold uppercase panel label (A, B, C ) to an axis."""
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=FONTSIZE['panel_label'], fontweight='bold',
            va='top', ha='left', color='#1A1A1A')


#  Statistical annotation 

def pvalue_to_stars(p):
    """Convert p-value to asterisk notation."""
    if p < 0.0001: return '****'
    if p < 0.001:  return '***'
    if p < 0.01:   return '**'
    if p < 0.05:   return '*'
    return 'ns'


def add_stat_bracket(ax, x1, x2, y_data_max, p, gap=0.04):
    """
    Draw a significance bracket above data between x-positions x1 and x2.
    y_data_max: current data maximum (used to position bracket).
    """
    ylim = ax.get_ylim()
    span = ylim[1] - ylim[0]
    y_line = y_data_max + span * gap
    y_text = y_line + span * 0.01

    ax.plot([x1, x1, x2, x2],
            [y_line - span * 0.01, y_line, y_line, y_line - span * 0.01],
            lw=0.8, color='#2C2C2C')

    label = pvalue_to_stars(p)
    fs = FONTSIZE['annotation'] if label != 'ns' else FONTSIZE['annotation'] - 0.5
    ax.text((x1 + x2) / 2, y_text, label,
            ha='center', va='bottom', fontsize=fs, color='#2C2C2C')

    # Expand y-axis to fit bracket
    ax.set_ylim(ylim[0], y_text + span * 0.08)


#  Violin + box + jitter (publication standard) 

def violin_with_stats(ax, data, x_col, y_col, palette=None, order=None,
                       point_size=2.0, alpha_violin=0.55, alpha_points=0.35):
    """
    Layered violin / boxplot / jitter with automatic Mann-Whitney U annotation.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    data : pandas.DataFrame
    x_col, y_col : str   column names
    palette : dict        {group: color}
    order : list          group display order
    """
    if palette is None:
        palette = get_condition_colors()
    if order is None:
        order = [k for k in palette if k in data[x_col].values]
    if not order:
        order = list(data[x_col].dropna().unique())

    # Violin
    vp = sns.violinplot(data=data, x=x_col, y=y_col, order=order,
                        palette=palette, inner=None, linewidth=0.75,
                        cut=0, ax=ax, saturation=0.9)
    for collection in ax.collections:
        if hasattr(collection, 'set_alpha'):
            collection.set_alpha(alpha_violin)

    # Narrow boxplot overlay
    sns.boxplot(data=data, x=x_col, y=y_col, order=order,
                width=0.10, palette=palette, linewidth=0.75,
                fliersize=0, ax=ax,
                boxprops=dict(alpha=0.9),
                medianprops=dict(color='white', linewidth=1.5),
                whiskerprops=dict(linewidth=0.75),
                capprops=dict(linewidth=0.75))

    # Jitter
    sns.stripplot(data=data, x=x_col, y=y_col, order=order,
                  palette=palette, size=point_size, alpha=alpha_points,
                  jitter=True, dodge=False, ax=ax, linewidth=0)

    # Statistical annotation
    if len(order) == 2:
        g1 = data[data[x_col] == order[0]][y_col].dropna().values
        g2 = data[data[x_col] == order[1]][y_col].dropna().values
        if len(g1) >= 3 and len(g2) >= 3:
            _, p = stats.mannwhitneyu(g1, g2, alternative='two-sided')
            y_max = data[y_col].max()
            add_stat_bracket(ax, 0, 1, y_max, p)

    sns.despine(ax=ax)
    return ax


#  UMAP axis formatting 

def format_umap_ax(ax, title=None):
    """Remove ticks from UMAP axis and optionally set title."""
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('UMAP 1', fontsize=FONTSIZE['axis_title'])
    ax.set_ylabel('UMAP 2', fontsize=FONTSIZE['axis_title'])
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    if title:
        ax.set_title(title, fontsize=FONTSIZE['title'], pad=4)


def add_umap_cell_type_labels(ax, adata, cell_type_col):
    """Add centered text labels at centroid of each cell type cluster."""
    umap = adata.obsm['X_umap']
    for ct in adata.obs[cell_type_col].unique():
        mask = adata.obs[cell_type_col] == ct
        cx, cy = umap[mask, 0].mean(), umap[mask, 1].mean()
        ax.text(cx, cy, ct, fontsize=5.5, ha='center', va='center',
                fontweight='bold', color='white',
                bbox=dict(boxstyle='round,pad=0.15', facecolor='#2C2C2C',
                          alpha=0.65, linewidth=0))


#  Volcano plot helper 

def volcano_plot(ax, result_df, label_col='names', fc_col='logfoldchanges',
                 pval_col='pvals_adj', fc_thresh=0.5, pval_thresh=0.05,
                 n_label=10, title=None):
    """
    Publication-quality volcano plot with adjustText gene labeling.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    result_df : pandas.DataFrame  scanpy DE output
    """
    from adjustText import adjust_text

    df = result_df.copy()
    df['log10p'] = -np.log10(df[pval_col].clip(lower=1e-300))

    # Color categories
    sig_up   = (df[pval_col] < pval_thresh) & (df[fc_col] >  fc_thresh)
    sig_down = (df[pval_col] < pval_thresh) & (df[fc_col] < -fc_thresh)
    ns       = ~(sig_up | sig_down)

    ax.scatter(df.loc[ns,   fc_col], df.loc[ns,   'log10p'],
               s=3, alpha=0.4, color='#CCCCCC', rasterized=True, linewidths=0)
    ax.scatter(df.loc[sig_up,   fc_col], df.loc[sig_up,   'log10p'],
               s=5, alpha=0.7, color=COLORS['LE'],     rasterized=True, linewidths=0)
    ax.scatter(df.loc[sig_down, fc_col], df.loc[sig_down, 'log10p'],
               s=5, alpha=0.7, color=COLORS['Normal'], rasterized=True, linewidths=0)

    # Threshold lines
    ax.axhline(-np.log10(pval_thresh), color='#888888', lw=0.6, ls='--')
    ax.axvline( fc_thresh,  color='#888888', lw=0.6, ls='--')
    ax.axvline(-fc_thresh,  color='#888888', lw=0.6, ls='--')

    # Label top genes
    top = pd.concat([
        df[sig_up].nlargest(n_label, fc_col),
        df[sig_down].nsmallest(n_label, fc_col)
    ])
    texts = []
    for _, row in top.iterrows():
        t = ax.text(row[fc_col], row['log10p'], row[label_col],
                    fontsize=FONTSIZE['annotation'], color='#1A1A1A')
        texts.append(t)
    if texts:
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle='-', color='#888888', lw=0.4))

    # Counts
    ax.text(0.98, 0.98, f'Up: {sig_up.sum():,}', transform=ax.transAxes,
            ha='right', va='top', color=COLORS['LE'], fontsize=FONTSIZE['annotation'])
    ax.text(0.02, 0.98, f'Down: {sig_down.sum():,}', transform=ax.transAxes,
            ha='left', va='top', color=COLORS['Normal'], fontsize=FONTSIZE['annotation'])

    ax.set_xlabel('Log fold change (LE / Normal)', fontsize=FONTSIZE['axis_title'])
    ax.set_ylabel('log(adjusted p-value)',        fontsize=FONTSIZE['axis_title'])
    if title:
        ax.set_title(title, fontsize=FONTSIZE['title'], pad=4)
    sns.despine(ax=ax)


#  Save 

def save_figure(fig, filepath, formats=('pdf', 'png'), dpi=300, tight=True):
    """
    Save in PDF (vector, journal submission) and PNG (preview).
    PDF uses fonttype 42 so text is editable in Adobe Illustrator.
    """
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)
    if tight:
        try:
            fig.tight_layout()
        except Exception:
            pass
    for fmt in formats:
        out = filepath.with_suffix(f'.{fmt}')
        fig.savefig(out, dpi=dpi, bbox_inches='tight',
                    format=fmt, facecolor='white', edgecolor='none')
        print(f'  Saved: {out.name}')


# Apply on import
import pandas as pd   # needed by volcano_plot
set_publication_style()


