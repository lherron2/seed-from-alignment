from setuptools import setup, find_packages

setup(
    name="sample_rna",          # e.g. "rnanneal"
    version="0.1.0",
    packages=find_packages(where="lib"),
    package_dir={"": "lib"},           # remove these two lines if not using a src/ layout
    python_requires=">=3.8",
    include_package_data=True,
)

