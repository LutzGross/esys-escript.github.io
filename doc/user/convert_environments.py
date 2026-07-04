#!/usr/bin/env python3
"""
Convert old Python documentation style LaTeX environments to standard LaTeX.

This script converts:
- \begin{classdesc}{Name}{args}...\end{classdesc}
- \begin{methoddesc}[Class]{method}{args}...\end{methoddesc}
- \begin{funcdesc}{func}{args}...\end{funcdesc}
- \begin{membdesc}[Class]{member}...\end{membdesc}
- \begin{datadesc}{name}...\end{datadesc}

To standard LaTeX that pandoc can understand.
"""

import re
import sys
from pathlib import Path


def convert_classdesc(content):
    r"""Convert classdesc environments.

    \begin{classdesc}{ClassName}{args} body \end{classdesc}
    ->
    \subsubsection*{class ClassName(args)}
    body
    """
    def replacer(match):
        classname = match.group(1)
        args = match.group(2).strip()
        body = match.group(3).strip()

        if args:
            sig = f"{classname}({args})"
        else:
            sig = classname

        return f"\\subsubsection*{{class {sig}}}\n{body}\n"

    pattern = r'\\begin\{classdesc\}\{([^}]+)\}\{([^}]*)\}(.*?)\\end\{classdesc\}'
    return re.sub(pattern, replacer, content, flags=re.DOTALL)


def convert_methoddesc(content):
    r"""Convert methoddesc environments.

    \begin{methoddesc}[Class]{method}{args} body \end{methoddesc}
    ->
    \paragraph{Class.method(args)}
    body
    """
    def replacer(match):
        classname = match.group(1) if match.group(1) else ''
        method = match.group(2)
        args = match.group(3).strip()
        body = match.group(4).strip()

        if classname:
            sig = f"{classname}.{method}({args})"
        else:
            sig = f"{method}({args})"

        return f"\\paragraph{{{sig}}}\n{body}\n"

    pattern = r'\\begin\{methoddesc\}(?:\[([^\]]*)\])?\{([^}]+)\}\{([^}]*)\}(.*?)\\end\{methoddesc\}'
    return re.sub(pattern, replacer, content, flags=re.DOTALL)


def convert_funcdesc(content):
    r"""Convert funcdesc environments.

    \begin{funcdesc}{func}{args} body \end{funcdesc}
    ->
    \paragraph{func(args)}
    body
    """
    def replacer(match):
        func = match.group(1)
        args = match.group(2).strip()
        body = match.group(3).strip()

        return f"\\paragraph{{{func}({args})}}\n{body}\n"

    pattern = r'\\begin\{funcdesc\}\{([^}]+)\}\{([^}]*)\}(.*?)\\end\{funcdesc\}'
    return re.sub(pattern, replacer, content, flags=re.DOTALL)


def convert_membdesc(content):
    r"""Convert membdesc environments.

    \begin{membdesc}[Class]{member} body \end{membdesc}
    ->
    \paragraph{Class.member}
    body
    """
    def replacer(match):
        classname = match.group(1) if match.group(1) else ''
        member = match.group(2)
        body = match.group(3).strip()

        if classname:
            sig = f"{classname}.{member}"
        else:
            sig = member

        return f"\\paragraph{{{sig}}}\n{body}\n"

    pattern = r'\\begin\{membdesc\}(?:\[([^\]]*)\])?\{([^}]+)\}(.*?)\\end\{membdesc\}'
    return re.sub(pattern, replacer, content, flags=re.DOTALL)


def convert_datadesc(content):
    r"""Convert datadesc environments.

    \begin{datadesc}{name} body \end{datadesc}
    ->
    \paragraph{name}
    body
    """
    def replacer(match):
        name = match.group(1)
        body = match.group(2).strip()

        return f"\\paragraph{{{name}}}\n{body}\n"

    pattern = r'\\begin\{datadesc\}\{([^}]+)\}(.*?)\\end\{datadesc\}'
    return re.sub(pattern, replacer, content, flags=re.DOTALL)


def convert_file(filepath, dry_run=False):
    """Convert a single LaTeX file."""
    with open(filepath, 'r', encoding='utf-8') as f:
        original = f.read()

    content = original

    # Count before
    counts_before = {
        'classdesc': len(re.findall(r'\\begin\{classdesc\}', content)),
        'methoddesc': len(re.findall(r'\\begin\{methoddesc\}', content)),
        'funcdesc': len(re.findall(r'\\begin\{funcdesc\}', content)),
        'membdesc': len(re.findall(r'\\begin\{membdesc\}', content)),
        'datadesc': len(re.findall(r'\\begin\{datadesc\}', content)),
    }

    total_before = sum(counts_before.values())
    if total_before == 0:
        return 0, {}

    # Apply conversions
    content = convert_classdesc(content)
    content = convert_methoddesc(content)
    content = convert_funcdesc(content)
    content = convert_membdesc(content)
    content = convert_datadesc(content)

    # Count after (should be 0)
    counts_after = {
        'classdesc': len(re.findall(r'\\begin\{classdesc\}', content)),
        'methoddesc': len(re.findall(r'\\begin\{methoddesc\}', content)),
        'funcdesc': len(re.findall(r'\\begin\{funcdesc\}', content)),
        'membdesc': len(re.findall(r'\\begin\{membdesc\}', content)),
        'datadesc': len(re.findall(r'\\begin\{datadesc\}', content)),
    }

    if not dry_run and content != original:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content)

    return total_before, counts_before


def main():
    import argparse

    parser = argparse.ArgumentParser(description='Convert Python doc style LaTeX to standard LaTeX')
    parser.add_argument('--dry-run', '-n', action='store_true',
                        help='Show what would be changed without modifying files')
    parser.add_argument('--file', '-f', type=str, default=None,
                        help='Convert a single file (default: all .tex files in doc/user/)')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Verbose output')
    args = parser.parse_args()

    if args.file:
        files = [Path(args.file)]
    else:
        # Find all .tex files in the same directory as this script
        script_dir = Path(__file__).parent
        files = list(script_dir.glob('*.tex'))

    total_converted = 0

    print(f"{'[DRY RUN] ' if args.dry_run else ''}Converting Python doc style environments to standard LaTeX\n")

    for filepath in sorted(files):
        count, details = convert_file(filepath, dry_run=args.dry_run)
        if count > 0:
            print(f"{filepath.name}: {count} environments converted")
            if args.verbose:
                for env, c in details.items():
                    if c > 0:
                        print(f"    {env}: {c}")
            total_converted += count

    print(f"\nTotal: {total_converted} environments {'would be ' if args.dry_run else ''}converted")

    if args.dry_run and total_converted > 0:
        print("\nRun without --dry-run to apply changes")


if __name__ == '__main__':
    main()
