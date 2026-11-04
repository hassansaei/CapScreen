import os
from pathlib import Path

from setuptools import setup, find_packages


def package_version():
    """Read the version stamped into capscreen/version.py."""
    version_file = Path(__file__).resolve().parent / "capscreen" / "version.py"
    for line in version_file.read_text().splitlines():
        if line.startswith("__version__"):
            return line.split("=", 1)[1].strip().strip('"').strip("'")
    raise RuntimeError(f"No __version__ assignment in {version_file}")


def install_version_hook():
    """Link the version bump hook in a local checkout."""
    root = Path(__file__).resolve().parent
    source = root / "githooks" / "post-commit"
    dest = root / ".git" / "hooks" / "post-commit"
    if not source.is_file() or not dest.parent.is_dir() or os.path.lexists(dest):
        return
    dest.symlink_to(Path("../../githooks/post-commit"))


install_version_hook()

setup(
    name="capscreen",
    version=package_version(),
    packages=find_packages(),
    package_data={
        "capscreen": [
            "config.json",
            "templates/*.html",
        ],
    },
    include_package_data=True,
    install_requires=[
        "setuptools>=67.0",
        "pathlib",
    ],
    entry_points={
        'console_scripts': [
            'capscreen=capscreen.cli:main',
        ],
    },
    author="Hassan Saei",
    author_email="hassan.saeiahan@gmail.com",
    description="A tool for screening AAV libraries for capsid variant detection in pull-down and transductions assays",
    python_requires=">=3.6",
)