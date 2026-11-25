"""Create release provenance; never build or publish an application binary."""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--engine-dir', type=Path, required=True)
    parser.add_argument('--sbom', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    def git(*arguments):
        return subprocess.check_output(['git', *arguments], cwd=root, text=True).strip()
    names = subprocess.check_output(['git', 'ls-files', '--cached', '--others', '--exclude-standard', '-z'], cwd=root).decode().split('\0')
    paths = [root / name for name in names if name]
    manifest = {
        'source_commit': git('rev-parse', 'HEAD'),
        'source_dirty': bool(git('status', '--porcelain')),
        'source_sha256': {str(path.relative_to(root)): digest(path) for path in paths
                          if path.is_file() and path.resolve() != args.output.resolve()},
        'engine_sha256': {path.name: digest(path) for path in sorted(args.engine_dir.iterdir()) if path.is_file()},
        'python_dependencies': {package: importlib.metadata.version(package)
                                for package in ('numpy', 'scipy', 'matplotlib', 'PyQt5', 'PyQt5-Qt5', 'PyQt5-sip')},
        'sbom_file': args.sbom.name,
        'sbom_sha256': digest(args.sbom),
        'note': 'A release must use a clean source commit. A dirty manifest is a local validation snapshot only.'
    }
    args.output.write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
