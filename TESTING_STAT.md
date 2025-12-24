# Testing Guide for stat.py with Group Roles Configuration

This guide explains how to test the statistical analysis module, including the new explicit group role configuration feature.

## Prerequisites

1. **Test Data Files:**
   - `merged.counts.tsv` - Merged count table with ID_WLG, peptide, and sample columns
   - `Sample_info.csv` - Sample metadata with group, tech_rep, bio_rep, batch_seq columns

2. **Current Test Data Groups:**
   Based on your Sample_info.csv, you have:
   - Input groups: `CD47_input`, `Pool_input`
   - Target groups: `CD47_WT`, `Pool_WT`
   - Background groups: `CD47_KO`, `Pool_KO`

## Testing Methods

### Method 1: Command Line Interface (CLI)

#### Test 1: Using Heuristics (No Config - Default Behavior)

This tests the fallback behavior when no explicit group roles are configured:

```bash
# From the project root directory
python -m capscreen.cli stat \
    --counts merged.counts.tsv \
    --sample-info Sample_info.csv \
    --output-dir test_output_heuristics \
    --verbose
```

**What to check:**
- Look for log messages like:
  - `"No explicit group role configuration found, using heuristics"`
  - `"Input groups: ['CD47_input', 'Pool_input']"`
  - `"Target groups: ['CD47_WT', 'Pool_WT']"`
  - `"Background groups: ['CD47_KO', 'Pool_KO']"`

#### Test 2: Using Explicit Group Roles Configuration

First, create a test config file with explicit group roles:

```bash
# Create a test config file
cat > test_config.json << 'EOF'
{
    "statistical_analysis": {
        "group_roles": {
            "input": ["CD47_input", "Pool_input"],
            "target": ["CD47_WT", "Pool_WT"],
            "background": ["CD47_KO", "Pool_KO"]
        },
        "analysis": {
            "n_cpus": 8,
            "lambda": 0.01,
            "low_expression_threshold": 10,
            "pca": {
                "n_components": 2,
                "random_state": 0,
                "top_n_features": 500
            }
        },
        "differential_expression": {
            "padj_threshold": 0.05,
            "log2fc_threshold": 1.0
        },
        "plots": {
            "scatter_highlight_top_n": 60,
            "bar_plot_top_n": 30,
            "dpi": {
                "standard": 300,
                "heatmap": 600
            },
            "figure_sizes": {
                "scatter": [6, 6],
                "bar": [6, 6],
                "rpm_scatter": [6, 4],
                "pca": [6, 5],
                "heatmap": [10, 10]
            },
            "scatter": {
                "all_points_size": 10,
                "all_points_alpha": 0.4,
                "highlighted_size": 40,
                "enable_highlighting": true
            },
            "pca": {
                "point_size": 120,
                "fontsize": 10
            },
            "fontsize": 10
        },
        "colors": {
            "target": "#d62728",
            "background_bead": "#1f77b4",
            "background_bead_fc": "#ff7f0e",
            "gray": "#808080",
            "highlight": "#d62728"
        }
    }
}
EOF
```

Then run with the config:

```bash
python -m capscreen.cli stat \
    --counts merged.counts.tsv \
    --sample-info Sample_info.csv \
    --output-dir test_output_explicit \
    --config test_config.json \
    --verbose
```

**What to check:**
- Look for log messages like:
  - `"Using explicit group role configuration from config file"`
  - `"Input groups: ['CD47_input', 'Pool_input']"`
  - `"Target groups: ['CD47_WT', 'Pool_WT']"`
  - `"Background groups: ['CD47_KO', 'Pool_KO']"`

#### Test 3: Testing with Misnamed Groups (Edge Case)

Create a config with atypically named groups to test that explicit config overrides heuristics:

```bash
# Create a config with groups that would be misclassified by heuristics
cat > test_config_misnamed.json << 'EOF'
{
    "statistical_analysis": {
        "group_roles": {
            "input": ["CD47_input", "Pool_input"],
            "target": ["CD47_WT", "Pool_WT"],
            "background": ["CD47_KO", "Pool_KO"]
        },
        "analysis": {
            "n_cpus": 8,
            "lambda": 0.01,
            "low_expression_threshold": 10,
            "pca": {
                "n_components": 2,
                "random_state": 0,
                "top_n_features": 500
            }
        }
    }
}
EOF
```

#### Test 4: Testing with Missing Groups (Error Handling)

Test what happens when config specifies groups that don't exist:

```bash
cat > test_config_missing.json << 'EOF'
{
    "statistical_analysis": {
        "group_roles": {
            "input": ["NonExistent_Input"],
            "target": ["CD47_WT"],
            "background": ["CD47_KO"]
        }
    }
}
EOF

python -m capscreen.cli stat \
    --counts merged.counts.tsv \
    --sample-info Sample_info.csv \
    --output-dir test_output_missing \
    --config test_config_missing.json \
    --verbose
```

**What to check:**
- Should see warning: `"Configured input groups not found in metadata: ['NonExistent_Input']"`
- Should see warning about unassigned groups being classified using heuristics

#### Test 5: Testing with Empty Config (Fallback to Heuristics)

Test that empty group_roles arrays fall back to heuristics:

```bash
cat > test_config_empty.json << 'EOF'
{
    "statistical_analysis": {
        "group_roles": {
            "input": [],
            "target": [],
            "background": []
        }
    }
}
EOF

python -m capscreen.cli stat \
    --counts merged.counts.tsv \
    --sample-info Sample_info.csv \
    --output-dir test_output_empty \
    --config test_config_empty.json \
    --verbose
```

**What to check:**
- Should see: `"No explicit group role configuration found, using heuristics"`

### Method 2: Direct Python Script Execution

You can also run stat.py directly:

```bash
python capscreen/scripts/stat.py \
    --counts merged.counts.tsv \
    --sample-info Sample_info.csv \
    --output-dir test_output_direct \
    --config test_config.json \
    --verbose
```

### Method 3: Programmatic Testing (Python)

Create a test script:

```python
# test_stat.py
from pathlib import Path
from capscreen.scripts import stat as stat_module

# Test with heuristics
print("Test 1: Using heuristics...")
success1 = stat_module.run_statistical_analysis(
    counts_file=Path("merged.counts.tsv"),
    sample_info_file=Path("Sample_info.csv"),
    output_dir=Path("test_output_prog1"),
    config_file=None,  # No config - use heuristics
    verbose=True
)
print(f"Result: {'Success' if success1 else 'Failed'}")

# Test with explicit config
print("\nTest 2: Using explicit config...")
success2 = stat_module.run_statistical_analysis(
    counts_file=Path("merged.counts.tsv"),
    sample_info_file=Path("Sample_info.csv"),
    output_dir=Path("test_output_prog2"),
    config_file=Path("test_config.json"),  # Use explicit config
    verbose=True
)
print(f"Result: {'Success' if success2 else 'Failed'}")
```

Run it:
```bash
python test_stat.py
```

## Expected Output Files

After successful execution, you should see these files in the output directory:

1. **Data Files:**
   - `normalized_counts.tsv`
   - `RPM_enrichment_analysis.tsv`
   - `correlation_statistics.tsv`
   - `DE_*.tsv` (differential expression results)

2. **Plots:**
   - `pca_all_features.png`
   - `pca_top500_features.png`
   - `pca_no_input_samples.png`
   - `correlation_scatterPlot.pdf`
   - `target_enrichment_scatter_*.png`
   - `top_target_specific_variants_*.png`
   - `rpm_scatter_*.png`

## Verification Checklist

- [ ] Logs show correct group classification (check for "Input groups:", "Target groups:", "Background groups:")
- [ ] When using config, logs show "Using explicit group role configuration from config file"
- [ ] When not using config, logs show "No explicit group role configuration found, using heuristics"
- [ ] All expected output files are created
- [ ] No errors in the log output
- [ ] RPM enrichment analysis includes correct groups
- [ ] Differential expression analysis runs for correct group pairs
- [ ] Plots are generated for target-background pairs

## Troubleshooting

1. **Import errors:** Make sure you're in the project root and dependencies are installed
2. **File not found:** Check that `merged.counts.tsv` and `Sample_info.csv` exist in the current directory
3. **Config not found:** If using `--config`, ensure the path is correct
4. **Group classification issues:** Check the logs to see which groups were identified and how

## Quick Test Command

Here's a quick one-liner to test with your existing data:

```bash
python -m capscreen.cli stat \
    --counts merged.counts.tsv \
    --sample-info Sample_info.csv \
    --output-dir test_output \
    --verbose
```

This will use heuristics and should complete successfully with your test data.

