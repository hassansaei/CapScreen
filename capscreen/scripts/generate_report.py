import os
import json
import argparse
from pathlib import Path
from jinja2 import Environment, FileSystemLoader, PackageLoader, ChoiceLoader
import re
import pandas as pd
from capscreen.version import __version__ as tool_version


def parse_bowtie2_log(filepath):
    if not os.path.exists(filepath):
        return {}
    with open(filepath) as f:
        lines = f.readlines()
    return {'raw': lines}


def parse_fastp_json(filepath):
    if not os.path.exists(filepath):
        return {}
    with open(filepath) as f:
        data = json.load(f)
    summary = data.get('summary', {})
    filtering_result = data.get('filtering_result', {})
    duplication = data.get('duplication', {})
    return {
        'before_filtering': summary.get('before_filtering', {}),
        'after_filtering': summary.get('after_filtering', {}),
        'filtering_result': filtering_result,
        'duplication': duplication
    }


def parse_pipeline_log(filepath):
    if not os.path.exists(filepath):
        return {}
    stats = {}
    merging_stats = {}
    variant_stats = {}
    with open(filepath) as f:
        lines = f.readlines()
    # Extract required lines
    for i, line in enumerate(lines):
        # Remove timestamp and log level for easier matching
        msg = re.sub(r'^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d+ - INFO - ', '', line)
        if 'Total reads in SAM file:' in msg:
            stats['total_reads'] = msg.split(':')[-1].strip().replace(',', '')
        elif 'Mapped reads:' in msg:
            stats['mapped_reads'] = msg.split(':')[-1].strip()
        elif 'Unmapped reads:' in msg:
            stats['unmapped_reads'] = msg.split(':')[-1].strip()
        elif 'Reads with flanking sequences:' in msg:
            stats['reads_with_flanking'] = msg.split(':')[-1].strip()
        elif 'Valid peptides after translation:' in msg:
            stats['valid_peptides'] = msg.split(':')[-1].strip()
        elif 'Unique peptides in reads:' in msg:
            merging_stats['unique_peptides_in_reads'] = msg.split(':')[-1].strip().replace(',', '')
        elif 'Unique peptides in reference:' in msg:
            merging_stats['unique_peptides_in_reference'] = msg.split(':')[-1].strip().replace(',', '')
        elif 'Peptides found in both:' in msg:
            merging_stats['peptides_found_in_both'] = msg.split(':')[-1].strip().replace(',', '')
        elif 'Peptides only in reads:' in msg:
            merging_stats['peptides_only_in_reads'] = msg.split(':')[-1].strip().replace(',', '')
        elif 'Total sequences:' in msg:
            variant_stats['total_sequences'] = msg.split(':')[-1].strip().replace(',', '')
        elif 'Assigned sequences:' in msg:
            variant_stats['assigned_sequences'] = msg.split(':')[-1].strip()
        elif 'Unassigned sequences:' in msg:
            variant_stats['unassigned_sequences'] = msg.split(':')[-1].strip()
        elif 'Unique variants found:' in msg:
            variant_stats['unique_variants_found'] = msg.split(':')[-1].strip().replace(',', '')
    return {
        'stats': stats,
        'merging_stats': merging_stats,
        'variant_stats': variant_stats
    }


def parse_pear_log(filepath):
    if not os.path.exists(filepath):
        return {}
    assembled = discarded = not_assembled = None
    with open(filepath) as f:
        for line in f:
            # Only match summary lines with counts and percentages
            if re.match(r'^Assembled reads\s+\.+: \d', line):
                assembled = line.split(':', 1)[-1].strip()
            elif re.match(r'^Discarded reads\s+\.+: \d', line):
                discarded = line.split(':', 1)[-1].strip()
            elif re.match(r'^Not assembled reads\s+\.+: \d', line):
                not_assembled = line.split(':', 1)[-1].strip()
    return {
        'assembled_reads': assembled,
        'discarded_reads': discarded,
        'not_assembled_reads': not_assembled
    }


def render_html_report(sample_name, fastp_data, pipeline_data, pear_data, pipeline_log_contents, output_path):
    template_dir = Path(__file__).parent.parent / 'templates'
    env = Environment(loader=FileSystemLoader(str(template_dir)))
    template = env.get_template('template_report.html')
    html = template.render(
        sample_name=sample_name,
        fastp=fastp_data,
        pipeline=pipeline_data,
        pear=pear_data,
        pipeline_log_contents=pipeline_log_contents
    )
    with open(output_path, 'w') as f:
        f.write(html)


def generate_report(sample_dir, output=None):
    sample_dir = Path(sample_dir)
    sample_name = sample_dir.name

    # Determine output filename
    if output is None:
        output_filename = f"{sample_name}_report.html"
    else:
        output_filename = output
    output_path = Path(output_filename)
    if not output_path.is_absolute():
        output_path = sample_dir / output_path

    # Read config
    config_path = sample_dir.parent.parent / 'config.json'
    if not config_path.exists():
        config_path = sample_dir.parent / 'config.json'
    if not config_path.exists():
        config_path = Path('capscreen/config.json')
    if config_path.exists():
        with open(config_path) as f:
            config = json.load(f)
    else:
        config = {}

    fastp_json = sample_dir / 'fastp.json'
    pipeline_log = sample_dir / f'{sample_name}.pipeline.log'
    pear_log = sample_dir / 'pear.pear.log'
    counts_tsv = sample_dir / f'{sample_name}.counts.tsv'

    fastp_data = parse_fastp_json(fastp_json)
    pipeline_data = parse_pipeline_log(pipeline_log)
    pear_data = parse_pear_log(pear_log)

    # Read full pipeline log for collapsible section
    if pipeline_log.exists():
        with open(pipeline_log) as f:
            pipeline_log_contents = f.read()
    else:
        pipeline_log_contents = ''

    # Prepare count matrix HTML and JSON for DataTables if enabled
    count_matrix_html = ''
    count_matrix_json = '[]'
    show_count_matrix = config.get('show_count_matrix', False)
    show_histogram = config.get('show_histogram', False)
    log2rpm_hist_path = None
    if show_count_matrix and counts_tsv.exists():
        df = pd.read_csv(counts_tsv, sep=None, engine='python')
        if 'ID_WLG' in df.columns:
            df = df.drop_duplicates(subset=['ID_WLG'])
            # Remove unwanted columns
            for col in ['insertions', 'deletions', 'matches', 'Unnamed: 0']:
                if col in df.columns:
                    df = df.drop(columns=[col])
            # Prepare HTML for first 10 rows
            count_matrix_html = df.head(10).to_html(index=False, escape=False, table_id='count-matrix-table')
            # Prepare JSON for all rows (for DataTables)
            count_matrix_json = df.to_json(orient='records')
            # Generate histogram for log2_RPM if present and show_histogram is True
            if show_histogram and 'log2_RPM' in df.columns:
                import matplotlib.pyplot as plt
                plt.figure(figsize=(6,4))
                plt.hist(df['log2_RPM'].dropna(), bins=30, color='#3498db', edgecolor='black')
                plt.xlabel('log2(RPM + 1)')
                plt.ylabel('Frequency')
                plt.title('Distribution of peptides log2(RPM + 1)')
                hist_path = sample_dir / 'log2rpm_hist.png'
                plt.tight_layout()
                plt.savefig(hist_path)
                plt.close()
                log2rpm_hist_path = hist_path.name

    # Resolve HTML template path for local, docker, and installed environments
    template_config = config.get('html_template')
    template_name = 'template_report.html'
    template_loaders = []
    searched_paths = []

    def try_template_path(candidate: Path) -> bool:
        nonlocal template_name
        if candidate is None:
            return False
        candidate = candidate.expanduser()
        if candidate.is_file():
            template_loaders.append(FileSystemLoader(str(candidate.parent)))
            template_name = candidate.name
            return True
        searched_paths.append(str(candidate))
        return False

    if template_config:
        cfg_path = Path(template_config)
        found = try_template_path(cfg_path)
        if not found and not cfg_path.is_absolute():
            repo_relative = (Path(__file__).parent.parent / cfg_path).resolve()
            found = try_template_path(repo_relative)
            if not found:
                cwd_relative = (Path.cwd() / cfg_path).resolve()
                found = try_template_path(cwd_relative)
        if not found:
            opt_path = Path("/opt/capscreen") / cfg_path
            try_template_path(opt_path)
    else:
        default_paths = [
            Path(__file__).parent.parent / 'templates' / template_name,
            Path.cwd() / 'capscreen' / 'templates' / template_name,
            Path('/opt/capscreen/capscreen/templates') / template_name,
        ]
        for path in default_paths:
            if try_template_path(path):
                break

    # Always fall back to package resources (works for installed distributions)
    template_loaders.append(PackageLoader('capscreen', 'templates'))
    env = Environment(loader=ChoiceLoader(template_loaders))

    try:
        template = env.get_template(template_name)
    except Exception as exc:
        raise FileNotFoundError(
            "Could not load HTML template. "
            f"Tried the following paths: {', '.join(searched_paths)}. "
            "Ensure the template exists or set 'html_template' to a valid file."
        ) from exc
    html = template.render(
        sample_name=sample_name,
        fastp=fastp_data,
        pipeline=pipeline_data,
        pear=pear_data,
        pipeline_log_contents=pipeline_log_contents,
        count_matrix_html=count_matrix_html,
        count_matrix_json=count_matrix_json,
        tool_version=tool_version,
        log2rpm_hist_path=log2rpm_hist_path
    )
    with open(output_path, 'w') as f:
        f.write(html)
    print(f'Report generated: {output_path}')


def main():
    parser = argparse.ArgumentParser(description='Generate HTML report for CapScreen pipeline.')
    parser.add_argument('sample_dir', help='Directory containing log and result files for the sample')
    parser.add_argument('--output', default=None, help='Output HTML file name')
    args = parser.parse_args()
    generate_report(args.sample_dir, args.output)


if __name__ == '__main__':
    main() 