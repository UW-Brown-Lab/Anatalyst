# sc_pipeline/modules/reportgenerator.py

import os
import logging
import markdown
import datetime
import yaml
from pathlib import Path
from sc_pipeline.core.module import AnalysisModule

class ReportGenerator(AnalysisModule):
    """Module for generating a markdown report from pipeline results."""

    # NOT ALL OF THESE ARE IMPLEMENTED YET but might be some day who knows
    PARAMETER_SCHEMA = {
        'report_file': {
            'type': str,
            'default': 'analysis_report.md',
            'description': 'Filename for the generated markdown report'
        },
        'generate_html': {
            'type': bool,
            'default': True,
            'description': 'Whether to also generate an HTML version of the report'
        },
        'include_pipeline_info': {
            'type': bool,
            'default': True,
            'description': 'Whether to include pipeline configuration information in the report'
        },
        'report_title': {
            'type': str,
            'default': None,
            'description': 'Custom title for the report. If not provided, uses the pipeline name'
        },
        'custom_css': {
            'type': str,
            'default': None,
            'description': 'Path to a custom CSS file to use for HTML styling'
        }
    }


    def __init__(self, name, params):
        super().__init__(name, params)
        self.logger = logging.getLogger(f"Module.{name}")
        self.required_inputs = []
        self.outputs = []

    def run(self, data_context):
        """
        Generate a Markdown report with all figures and analysis results.
        
        Args:
            data_context: The shared data context
            
        Returns:
            bool: True if successful, False otherwise
        """
        try:
            # Get output directory
            output_dir = data_context.get("OUTPUT_DIR", os.getcwd())
            report_file = self.params.get("report_file", "analysis_report.md")
            report_path = os.path.join(output_dir, report_file)
            
            # Create report's directory if it doesn't exist
            reports_dir = os.path.dirname(report_path)
            os.makedirs(reports_dir, exist_ok=True)
            
            # Get report figures data
            report_figures = data_context.get("REPORT_FIGURES", {})
            
            # Get pipeline metadata if available
            pipeline_config = data_context.get("CONFIG", None)
            
            # Start building the report
            self.logger.info(f"Generating report at: {report_path}")
            
            # Create the markdown content
            md_content = self._generate_markdown(pipeline_config, report_figures, output_dir)
            
            # Write the markdown file
            with open(report_path, 'w') as f:
                f.write(md_content)
            
            # Generate HTML report if requested
            if self.params.get("generate_html", True):
                html_path = os.path.splitext(report_path)[0] + ".html"
                self._generate_html(md_content, html_path)
            
            self.logger.info(f"Report generation completed successfully")
            return True
            
        except Exception as e:
            self.logger.error(f"Error generating report: {e}", exc_info=True)
            return False
        
    def _generate_markdown(self, pipeline_config, report_figures, output_dir):
        """Generate the markdown content for the report."""
        
        # Start with a header
        now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        md = [
            f"# {pipeline_config['pipeline']['name']} Report",
            f"*Generated on: {now}*",
            "",
            "## Overview",
            "",
            "This report summarizes the results of the single-cell RNA-seq analysis pipeline.",
            "",
            "## Analysis Results",
            ""
        ]
        

        # Add figures organized by module
        if not report_figures:
            md.append("No analysis figures were generated during the pipeline run.")
        else:
            # Process each module's figures
            for module_name, figures in report_figures.items():
                # Add a section for this module
                md.append(f"### {module_name}")
                md.append("")
                
                # Add each figure with its metadata
                for i, fig in enumerate(figures):
                    # Add figure title if available
                    if fig.get('title'):
                        md.append(f"#### {fig['title']}")
                        md.append("")
                    
                    # Add description if available
                    if fig.get('description'):
                        md.append(f"{fig['description']}")
                        md.append("")
                    
                    # Add the image if available
                    if fig.get('image_path'):
                        # Create a relative path for the image
                        img_path = fig['image_path']
                        rel_path = os.path.relpath(img_path, output_dir)
                        md.append(f"![{fig.get('title', f'Figure {i+1}')}]({rel_path})")
                        md.append("")
                    
                    # Add caption if available
                    if fig.get('caption'):
                        md.append(f"*{fig['caption']}*")
                        md.append("")
                
                # Add a separator between modules
                md.append("---")
                md.append("")
        

        md.append("## Pipeline Configuration")
        md.append("")
        md.append("The following configuration was used for this pipeline run:")
        md.append("")
        md.append("```yaml")
        try:
            config_yaml = yaml.dump(pipeline_config, default_flow_style=False, sort_keys=False)
            md.append(config_yaml)
        except Exception as e:
            self.logger.error(f"Error dumping config to YAML: {e}")
            md.append("# Error: Could not generate YAML configuration")
            
        md.append("```")
        md.append("")

        return "\n".join(md)
    
    def _generate_html(self, md_content, html_path):
        """Generate an HTML report from markdown content."""
        try:
            # Convert markdown to HTML
            html = markdown.markdown(md_content, extensions=[
                'markdown.extensions.tables',
                'markdown.extensions.fenced_code',
                'markdown.extensions.codehilite',
                'markdown.extensions.toc'
            ])
            
            # Add basic styling
            styled_html = f"""
            <!DOCTYPE html>
            <html>
            <head>
                <meta charset="UTF-8">
                <meta name="viewport" content="width=device-width, initial-scale=1.0">
                <title>scRNA-seq Analysis Report</title>
                <style>
                    body {{
                        font-family: Arial, sans-serif;
                        line-height: 1.6;
                        color: #333;
                        max-width: 1200px;
                        margin: 0 auto;
                        padding: 20px;
                    }}
                    h1, h2, h3, h4 {{
                        color: #2c3e50;
                    }}
                    img {{
                        max-width: 100%;
                        height: auto;
                        border: 1px solid #ddd;
                        border-radius: 4px;
                        padding: 5px;
                    }}
                    code {{
                        background-color: #f5f5f5;
                        padding: 2px 4px;
                        border-radius: 4px;
                    }}
                    pre {{
                        background-color: #f5f5f5;
                        padding: 10px;
                        border-radius: 4px;
                        overflow-x: auto;
                    }}
                    table {{
                        border-collapse: collapse;
                        width: 100%;
                    }}
                    th, td {{
                        border: 1px solid #ddd;
                        padding: 8px;
                    }}
                    tr:nth-child(even) {{
                        background-color: #f2f2f2;
                    }}
                    th {{
                        padding-top: 12px;
                        padding-bottom: 12px;
                        text-align: left;
                        background-color: #2c3e50;
                        color: white;
                    }}
                </style>
            </head>
            <body>
                {html}
            </body>
            </html>
            """
            
            # Write the HTML file
            with open(html_path, 'w') as f:
                f.write(styled_html)
                
            return True
            
        except Exception as e:
            self.logger.error(f"Error generating HTML report: {e}", exc_info=True)
            return False