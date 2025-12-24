#!/bin/bash
# Simple test script for stat.py
# This script tests both heuristic and explicit config-based group identification

set -e  # Exit on error

echo "=========================================="
echo "Testing stat.py - Group Roles Configuration"
echo "=========================================="
echo ""

# Check if required files exist
if [ ! -f "merged.counts.tsv" ]; then
    echo "ERROR: merged.counts.tsv not found!"
    exit 1
fi

if [ ! -f "Sample_info.csv" ]; then
    echo "ERROR: Sample_info.csv not found!"
    exit 1
fi

# Test 1: Heuristics (no config)
echo "Test 1: Testing with heuristics (no config)..."
echo "----------------------------------------"
python -m capscreen.cli stat \
    --counts merged.counts.tsv \
    --sample-info Sample_info.csv \
    --output-dir test_output_heuristics \
    --verbose 2>&1 | grep -E "(Input groups|Target groups|Background groups|explicit|heuristics)" || true
echo "✓ Test 1 completed"
echo ""

# Test 2: Explicit config
echo "Test 2: Testing with explicit group_roles config..."
echo "----------------------------------------"

# Create test config
cat > test_config_temp.json << 'EOF'
{
    "statistical_analysis": {
        "group_roles": {
            "input": ["CD47_input", "Pool_input"],
            "target": ["CD47_WT", "Pool_WT"],
            "background": ["CD47_KO", "Pool_KO"]
        },
        "analysis": {
            "n_cpus": 4,
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
            "dpi": {"standard": 300, "heatmap": 600},
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
            "pca": {"point_size": 120, "fontsize": 10},
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

python -m capscreen.cli stat \
    --counts merged.counts.tsv \
    --sample-info Sample_info.csv \
    --output-dir test_output_explicit \
    --config test_config_temp.json \
    --verbose 2>&1 | grep -E "(Input groups|Target groups|Background groups|explicit|heuristics)" || true
echo "✓ Test 2 completed"
echo ""

# Cleanup
rm -f test_config_temp.json

echo "=========================================="
echo "Testing complete!"
echo "=========================================="
echo ""
echo "Check output directories:"
echo "  - test_output_heuristics/"
echo "  - test_output_explicit/"
echo ""
echo "Both should contain the same analysis results, but may show different"
echo "group identification methods in the logs."

