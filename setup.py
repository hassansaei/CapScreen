from setuptools import setup, find_packages

setup(
    name="capscreen",
    version="1.0.0-alpha.6",
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