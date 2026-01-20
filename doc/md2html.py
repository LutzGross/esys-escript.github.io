#!/usr/bin/env python3
"""
Convert markdown files to HTML with Sphinx theme
"""

import sys
import os

try:
    import markdown
except ImportError:
    print("ERROR: markdown library not found. Please install python3-markdown")
    sys.exit(1)

def get_sphinx_template(title, content, current_page=""):
    """Generate HTML using Sphinx classic theme structure"""

    # Determine navigation links based on current page
    nav_links = []
    if current_page != "index":
        nav_links.append('<li class="right" style="margin-right: 10px"><a href="index.html">Home</a> |</li>')
    if current_page != "installation":
        nav_links.append('<li class="right"><a href="installation.html">Installation</a> |</li>')
    nav_links.append('<li class="right"><a href="sphinx_api/index.html">Python API</a> |</li>')
    nav_links.append('<li class="right"><a href="user/user.pdf">User Guide (PDF)</a></li>')

    nav_html = '\n        '.join(nav_links)

    template = """<!DOCTYPE html>
<html lang="en" data-content_root="./">
  <head>
    <meta charset="utf-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1.0" />
    <title>{title} &#8212; esys-escript documentation</title>
    <link rel="stylesheet" type="text/css" href="sphinx_api/_static/pygments.css" />
    <link rel="stylesheet" type="text/css" href="sphinx_api/_static/classic.css" />
    <link rel="stylesheet" type="text/css" href="sphinx_api/_static/custom.css" />
  </head>
  <body>
    <div class="related" role="navigation" aria-label="Related">
      <h3>Navigation</h3>
      <ul>
        {nav_links}
        <li class="nav-item nav-item-0"><a href="index.html">esys-escript documentation</a> &#187;</li>
      </ul>
    </div>

    <div class="document">
      <div class="documentwrapper">
        <div class="bodywrapper">
          <div class="body" role="main">
            {content}
            <div class="clearer"></div>
          </div>
        </div>
      </div>
      <div class="sphinxsidebar" role="navigation" aria-label="Main">
        <div class="sphinxsidebarwrapper">
          <h3><a href="index.html">Table of Contents</a></h3>
          <ul>
            <li><a href="index.html">Documentation Home</a></li>
            <li><a href="installation.html">Installation Guide</a></li>
            <li><a href="sphinx_api/index.html">Python API Reference</a></li>
            <li><a href="user/user.pdf">User Guide (PDF)</a></li>
          </ul>
          <h3>Downloads</h3>
          <ul>
            <li><a href="escript_examples.zip">Examples (ZIP)</a></li>
            <li><a href="escript_examples.tar.gz">Examples (TAR.GZ)</a></li>
          </ul>
          <h3>Resources</h3>
          <ul>
            <li><a href="https://github.com/esys-escript/esys-escript.github.io">GitHub Repository</a></li>
            <li><a href="https://github.com/esys-escript/esys-escript.github.io/issues">Report Issues</a></li>
          </ul>
        </div>
      </div>
    </div>

    <div class="footer" role="contentinfo">
      <p>&copy; Copyright 2003-2026, esys.escript group.</p>
      <p>
        Last updated on {date}.
        Created using <a href="https://www.python.org/">Python</a> and
        <a href="https://www.sphinx-doc.org/">Sphinx</a> theme.
      </p>
    </div>
  </body>
</html>
"""

    from datetime import datetime
    date = datetime.now().strftime("%b %d, %Y")

    return template.format(
        title=title,
        content=content,
        nav_links=nav_html,
        date=date
    )

def convert_readme_to_html(md_file, html_file):
    """Convert options.md to index.html with Sphinx theme"""

    with open(md_file, 'r') as f:
        md_content = f.read()

    # Convert markdown to HTML
    html_body = markdown.markdown(md_content, extensions=['extra', 'tables'])

    # Fix links to point to the correct HTML files
    html_body = html_body.replace('href="installation.md"', 'href="installation.html"')
    html_body = html_body.replace('href="./user.pdf"', 'href="user/user.pdf"')

    # Wrap in section tags for Sphinx styling
    html_body = '<section id="esys-escript">\n' + html_body + '\n</section>'

    # Create the documentation resources navigation section
    doc_nav_section = """<section id="documentation-resources">
<h2>Documentation Resources</h2>
<p><strong>Installation Guide</strong></p>
<ul class="simple">
<li><a class="reference external" href="installation.html">Installation Guide (HTML)</a></li>
<li>Source: <code>installation.md</code> in the repository root</li>
</ul>
<p><strong>User Guide (PDF)</strong></p>
<ul class="simple">
<li><a class="reference external" href="user/user.pdf">User Guide PDF</a></li>
</ul>
<p><strong>Python API Reference</strong></p>
<ul class="simple">
<li><a class="reference external" href="sphinx_api/index.html">Python API Documentation</a></li>
</ul>
<p><strong>Examples</strong></p>
<ul class="simple">
<li><a class="reference external" href="escript_examples.zip">Example Scripts (ZIP)</a></li>
<li><a class="reference external" href="escript_examples.tar.gz">Example Scripts (TAR.GZ)</a></li>
</ul>
</section>
"""

    # Insert documentation resources section
    if '<h2>Funding' in html_body or '<h2>The project was funded by the</h2>' in html_body:
        html_body = html_body.replace('<h2>Funding', doc_nav_section + '<h2>Funding')
        html_body = html_body.replace('<h2>The project was funded by the</h2>',
                                     doc_nav_section + '<h2>The project was funded by the</h2>')
    else:
        # Fallback: add before closing section tag
        html_body = html_body.replace('</section>', doc_nav_section + '</section>')

    full_html = get_sphinx_template("esys-escript", html_body, "index")

    with open(html_file, 'w') as f:
        f.write(full_html)

    print(f"Converted {md_file} to {html_file}")

def convert_installation_to_html(md_file, html_file):
    """Convert installation.md to HTML with Sphinx theme"""

    with open(md_file, 'r') as f:
        md_content = f.read()

    # Convert markdown to HTML
    html_body = markdown.markdown(md_content, extensions=['extra', 'tables', 'fenced_code', 'toc'])

    # Fix links to point to the correct HTML files
    html_body = html_body.replace('href="./scons/templates/options.md"', 'href="options.html"')

    # Wrap in section tag for Sphinx styling
    html_body = '<section id="installation-guide">\n' + html_body + '\n</section>'

    full_html = get_sphinx_template("Installation Guide", html_body, "installation")

    with open(html_file, 'w') as f:
        f.write(full_html)

    print(f"Converted {md_file} to {html_file}")

def convert_options_to_html(md_file, html_file):
    """Convert options.md to HTML with Sphinx theme"""

    with open(md_file, 'r') as f:
        md_content = f.read()

    # Convert markdown to HTML
    html_body = markdown.markdown(md_content, extensions=['extra', 'tables'])

    full_html = get_sphinx_template("Build Options Reference", html_body, "options")

    with open(html_file, 'w') as f:
        f.write(full_html)

    print(f"Converted {md_file} to {html_file}")

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Usage: md2html.py <type> <input.md> <output.html>")
        print("  type: 'readme', 'installation', or 'options'")
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
    elif doc_type == 'options':
        convert_options_to_html(input_file, output_file)
    else:
        print(f"ERROR: Unknown type '{doc_type}'. Use 'readme', 'installation', or 'options'")
        sys.exit(1)
