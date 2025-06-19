from setuptools import setup, find_packages

setup(
    name="capscreen",
    version="1.0.0-alpha.2",
    packages=find_packages(),
    package_data={
        "capscreen": [
            "config.json",
        ],
    },
    include_package_data=True,
    install_requires=[
        "pathlib",
    ],
    entry_points={
        'console_scripts': [
            'capscreen=capscreen.CapScreen:main',
        ],
    },
    author="Hassan Saei",
    author_email="hassan.saeiahan@gmail.com",
    description="A tool for screening AAV libraries for capsid variant detection in pull-down and transductions assays",
    python_requires=">=3.6",
)