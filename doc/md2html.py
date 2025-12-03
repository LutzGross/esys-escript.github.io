#!/usr/bin/env python3
"""
Convert markdown files to HTML with custom templates
"""

import sys
import os

try:
    import markdown
except ImportError:
    print("ERROR: markdown library not found. Please install python3-markdown")
    sys.exit(1)

def convert_readme_to_html(md_file, html_file):
    """Convert README.md to index.html with links to sphinx and installation docs"""

    with open(md_file, 'r') as f:
        md_content = f.read()

    # Convert markdown to HTML
    html_body = markdown.markdown(md_content, extensions=['extra', 'tables'])

    # Fix links to point to the correct HTML files
    html_body = html_body.replace('href="installation.md"', 'href="installation.html"')
    html_body = html_body.replace('href="./user.pdf"', 'href="user/user.pdf"')

    # Create the documentation resources navigation section
    doc_nav_section = """<h2>Documentation Resources</h2>
<ul>
<li><a href="sphinx_api/index.html">Python API Reference (Sphinx)</a></li>
<li><a href="installation.html">Installation Guide</a></li>
<li><a href="user/user.pdf">User Guide (PDF)</a></li>
<li><a href="escript_examples.zip">Example Scripts (ZIP)</a></li>
<li><a href="escript_examples.tar.gz">Example Scripts (TAR.GZ)</a></li>
</ul>
"""

    # Insert documentation resources section after "Using esys-escript" section
    # Look for the heading that comes after "Using esys-escript"
    if '<h2>The project was funded by the</h2>' in html_body:
        html_body = html_body.replace('<h2>The project was funded by the</h2>',
                                     doc_nav_section + '<h2>The project was funded by the</h2>')
    else:
        # Fallback: add it at the end if we can't find the marker
        html_body = html_body + doc_nav_section

    # Create full HTML page with navigation
    html_template = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>esys-escript Documentation</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
            line-height: 1.6;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            color: #333;
        }}
        h1 {{
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            margin-top: 30px;
            border-bottom: 1px solid #bdc3c7;
            padding-bottom: 5px;
        }}
        a {{
            color: #3498db;
            text-decoration: none;
        }}
        a:hover {{
            text-decoration: underline;
        }}
        code {{
            background: #f4f4f4;
            padding: 2px 6px;
            border-radius: 3px;
            font-family: "Courier New", Courier, monospace;
        }}
        pre {{
            background: #f4f4f4;
            padding: 15px;
            border-radius: 5px;
            overflow-x: auto;
        }}
        .doc-nav {{
            background: #ecf0f1;
            padding: 20px;
            border-radius: 5px;
            margin: 30px 0;
        }}
        .doc-nav h3 {{
            margin-top: 0;
            color: #2c3e50;
        }}
        .doc-nav ul {{
            list-style: none;
            padding: 0;
        }}
        .doc-nav li {{
            margin: 10px 0;
            font-size: 1.1em;
        }}
        .doc-nav li:before {{
            content: "📄 ";
            margin-right: 8px;
        }}
    </style>
</head>
<body>
    {content}
</body>
</html>
"""

    full_html = html_template.format(content=html_body)

    with open(html_file, 'w') as f:
        f.write(full_html)

    print(f"Converted {md_file} to {html_file}")

def convert_installation_to_html(md_file, html_file):
    """Convert installation.md to HTML with navigation"""

    with open(md_file, 'r') as f:
        md_content = f.read()

    # Convert markdown to HTML
    html_body = markdown.markdown(md_content, extensions=['extra', 'tables', 'fenced_code'])

    # Create full HTML page with navigation
    html_template = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>esys-escript Installation Guide</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
            line-height: 1.6;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            color: #333;
        }}
        h1 {{
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            margin-top: 30px;
            border-bottom: 1px solid #bdc3c7;
            padding-bottom: 5px;
        }}
        h3 {{
            color: #34495e;
            margin-top: 20px;
        }}
        a {{
            color: #3498db;
            text-decoration: none;
        }}
        a:hover {{
            text-decoration: underline;
        }}
        code {{
            background: #f4f4f4;
            padding: 2px 6px;
            border-radius: 3px;
            font-family: "Courier New", Courier, monospace;
        }}
        pre {{
            background: #2c3e50;
            color: #ecf0f1;
            padding: 15px;
            border-radius: 5px;
            overflow-x: auto;
        }}
        pre code {{
            background: transparent;
            color: inherit;
            padding: 0;
        }}
        ul {{
            margin: 10px 0;
        }}
        li {{
            margin: 5px 0;
        }}
        .back-link {{
            background: #ecf0f1;
            padding: 10px;
            border-radius: 5px;
            margin-bottom: 20px;
            display: inline-block;
        }}
    </style>
</head>
<body>
    <div class="back-link">
        <a href="index.html">← Back to Documentation Home</a>
    </div>

    {content}

    <div class="back-link" style="margin-top: 30px;">
        <a href="index.html">← Back to Documentation Home</a>
    </div>
</body>
</html>
"""

    full_html = html_template.format(content=html_body)

    with open(html_file, 'w') as f:
        f.write(full_html)

    print(f"Converted {md_file} to {html_file}")

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Usage: md2html.py <type> <input.md> <output.html>")
        print("  type: 'readme' or 'installation'")
        sys.exit(1)

    doc_type = sys.argv[1]
    input_file = sys.argv[2]
    output_file = sys.argv[3]

    if not os.path.exists(input_file):
        print(f"ERROR: Input file {input_file} does not exist")
        sys.exit(1)

    # Ensure output directory exists
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if doc_type == 'readme':
        convert_readme_to_html(input_file, output_file)
    elif doc_type == 'installation':
        convert_installation_to_html(input_file, output_file)
    else:
        print(f"ERROR: Unknown type '{doc_type}'. Use 'readme' or 'installation'")
        sys.exit(1)